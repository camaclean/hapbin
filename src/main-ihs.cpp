/*
 * Hapbin: A fast binary implementation EHH, iHS, and XPEHH
 * Copyright (C) 2014  Colin MacLean <s0838159@sms.ed.ac.uk>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "hapbin.hpp"
#include "hapmap.hpp"
#include "ihsfinder.hpp"
#include "calcselect.hpp"
#include "config.h"
#include "argparse.hpp"

#if MPI_FOUND
#include <mpi.h>
#endif
#include <functional>
#include <cstdlib>

#ifdef HAVE_THRUST
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/for_each.h>
#include <thrust/sequence.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/scan.h>
#include <thrust/set_operations.h>
#endif

struct ParentMetadata
{
    std::size_t snp;
    std::ptrdiff_t offset;
    std::size_t parent_index;
    bool direction;
};

struct ChildMetadata
{
    std::size_t snp;
    std::ptrdiff_t offset;
    std::size_t child_index;
    unsigned long long binom;
    bool direction;
    bool hit_end;
};

struct SnpMetadata
{
    std::size_t snp;
    unsigned long long sum;
    unsigned long long core_binom;
    double ihh0;
    double ihh1;
};

int main(int argc, char** argv)
{
    int ret = 0;
    std::size_t numSnps;
#if MPI_FOUND
    MPI_Init(&argc, &argv);
#endif
    Argument<bool> help('h', "help", "Show this help", true, false);
    Argument<bool> version('v', "version", "Version information", true, false);
    Argument<std::string> hap('d', "hap", "Hap file", false, false, "");
    Argument<std::string> map('m', "map", "Map file", false, false, "");
    Argument<double> cutoff('c', "cutoff", "EHH cutoff value (default: 0.05)", false, false, 0.05);
    Argument<double> minMAF('f', "minmaf", "Minimum allele frequency (default: 0.05)", false, false, 0.05);
    Argument<int> binfac('b', "bin", "Number of frequency bins for iHS normalization (default: 50)", false, false, 50);
    Argument<unsigned long long> scale('s', "scale", "Gap scale parameter in bp, used to scale gaps > scale parameter as in Voight, et al.", false, false, 20000);
    Argument<bool> binom('a', "binom", "Use binomial coefficients rather than frequency squared for EHH", true, false);
    Argument<unsigned long long> maxExtend('e', "max-extend", "Maximum distance in bp to traverse when calculating EHH (default: 0 (disabled))", false, false, 0);
    Argument<std::string> outfile('o', "out", "Output file", false, false, "out.txt");
    ArgParse argparse({&help, &version, &hap, &map, &outfile, &cutoff, &minMAF, &scale, &binfac, &maxExtend, &binom}, "Usage: ihsbin --map input.map --hap input.hap [--ascii] [--out outfile]");
    if (!argparse.parseArguments(argc, argv)) 
    {
        ret = 1;
        goto out;
    }
    if (help.value())
    {
        argparse.showHelp();
        goto out;
    }
    else if (version.value())
    {
        argparse.showVersion();
        goto out;
    }
    else if (!hap.wasFound() || !map.wasFound())
    {
        std::cout << "Please specify --hap and --map." << std::endl;
        ret = 2;
        goto out;
    }
    
    numSnps = HapMap::querySnpLength(hap.value().c_str());
    std::cout << "Chromosomes per SNP: " << numSnps << std::endl;
 out:
    HapMap m;
    m.loadHap(hap.value().c_str());
    m.loadMap(map.value().c_str());
    std::size_t num_snps    = m.numSnps();
    std::size_t snpDataSize = m.snpDataSize();
    std::size_t snpDataSizeULL = m.snpDataSizeULL();
    unsigned long long mask = ::bitsetMask<unsigned long long>(snpDataSizeULL);
    std::size_t max_breadth = 200000;
    std::size_t start_size = 200;
    std::size_t num_leaves = start_size;
    std::size_t active_snps_prev = start_size;
    typedef typename thrust::device_vector<unsigned long long>::iterator ULLIterator;
    typedef typename thrust::tuple<ULLIterator, ULLIterator, ULLIterator> ULLIteratorTuple;
    typedef typename thrust::zip_iterator<ULLIteratorTuple> ULL3Iterator;
    typedef typename thrust::device_vector<std::size_t>::iterator SizetIterator;
    typedef typename thrust::device_vector<std::ptrdiff_t>::iterator PtrdifftIterator;
    typedef typename thrust::tuple<SizetIterator, PtrdifftIterator, SizetIterator, SizetIterator> SizetPtrdifftTuple;
    using SizetPtrdifftIterator = thrust::zip_iterator<SizetPtrdifftTuple>;

    unsigned long long* hd = (unsigned long long*) m.rawData();
    thrust::device_vector<unsigned long long> map_data(m.numSnps()*snpDataSizeULL);
    for (std::size_t i = 0; i < num_snps; ++i)
        thrust::copy_n(&hd[i*snpDataSize*VEC],snpDataSizeULL, map_data.begin()+i*snpDataSizeULL);
    thrust::device_vector<double>             genetic_distance(num_snps);
    thrust::device_vector<unsigned long long> physical_distance(num_snps);
    thrust::copy_n(m.geneticPositionData(), num_snps, genetic_distance.begin());
    thrust::copy_n(m.physicalPositionData(), num_snps, physical_distance.begin());

    thrust::device_vector<unsigned long long> parent0(snpDataSizeULL*(max_breadth+1));
    thrust::device_vector<unsigned long long> parent1(snpDataSizeULL*(max_breadth+1));
    thrust::device_vector<unsigned long long> child0(snpDataSizeULL*(2*max_breadth));
    thrust::device_vector<unsigned long long> child1(snpDataSizeULL*(2*max_breadth));
    thrust::device_vector<unsigned long long> sums(num_snps);
    thrust::device_vector<unsigned long long> binom_core(num_snps);

    thrust::device_vector<ParentMetadata> parent_metadata(max_breadth);
    thrust::device_vector<ChildMetadata>  child_metadata(2*max_breadth);
    thrust::device_vector<std::size_t>    snps(max_breadth);
    thrust::device_vector<std::size_t>    reduced_snps(max_breadth);
    thrust::device_vector<std::size_t>    reduced_snps_prev(max_breadth);
    thrust::device_vector<std::size_t>    reduced_snps_prev_intersection(max_breadth);
    thrust::device_vector<std::size_t>    snp_branches(2*max_breadth);
    thrust::device_vector<unsigned long long> binoms(max_breadth);
    thrust::device_vector<unsigned long long> reduced_binoms(max_breadth);
    thrust::device_vector<float>          ehh0_prev(max_breadth);
    thrust::device_vector<float>          ehh0_prev_intersection(max_breadth);
    thrust::device_vector<float>          ehh1_prev(max_breadth);
    thrust::device_vector<float>          ehh1_prev_reduced(max_breadth);
    thrust::device_vector<float>          ehh0(num_snps);
    thrust::device_vector<float>          ehh1(num_snps);
    thrust::device_vector<float>         maf(num_snps);
    //
    //thrust::device_vector<std::ptrdiff_t> offsets(max_breadth+1);
    //thrust::device_vector<std::ptrdiff_t> offset_branches(2*max_breadth+1);
    //thrust::device_vector<std::size_t>    parent_indexes(max_breadth+1);
    //thrust::device_vector<std::size_t>    branch_indexes(2*max_breadth+1);
    //thrust::device_vector<unsigned int>   branch_counts(2*max_breadth+1);
    thrust::device_vector<std::size_t>    nonzero_branches(2*max_breadth);
    thrust::device_vector<std::size_t>    compact_branches(2*max_breadth);
    //thrust::device_vector<bool>           direction(max_breadth+1);
    //thrust::device_vector<bool>           hit_end(max_breadth+1);

    //Set parent0 and parent1 to be
    //parent0:
    //{junk,core0_0,core0_1,core0_2,core0_3,...}
    //parent1:
    //{junk,core1_0,core1_1,core1_2,core1_3,...}
    auto parent0it = parent0.begin();
    auto parent1it = parent1.begin();
    for (std::size_t i = 0; i < start_size; ++i)
    {
        parent0it += snpDataSizeULL;
        parent1it += snpDataSizeULL;
        thrust::copy_n(&hd[i*snpDataSizeULL], snpDataSizeULL, parent0it);
        thrust::copy_n(&hd[i*snpDataSizeULL], snpDataSizeULL, parent1it);
    }

    thrust::for_each_n(thrust::make_counting_iterator<std::size_t>(0), num_snps,
                       [&](std::size_t i){
        int count = 0;
        for(std::size_t j = 0; j < snpDataSizeULL; ++j)
            count += popcount1(map_data[i*snpDataSizeULL+j]);
        maf[i] = 1.0/(float)count;
        binom_core[i] = count*(count-1)/2;
    });

    thrust::for_each_n(thrust::make_counting_iterator<std::size_t>(0), start_size,
                       [&](std::size_t i){
        ParentMetadata m;
        m.snp = i;
        m.offset = 1;
        m.parent_index = i+1;
        parent_metadata[i] = m;
        reduced_snps_prev[i] = i;
        ehh0_prev[i] = 1.0;
        ehh1_prev[i] = 1.0;
    });
    //Set snp indexes to be:
    //{0, 1, 2, 3, 4, ..., start_size, __SIZE_MAX___, ...}
    //thrust::fill(snps.begin()+start_size, snps.end(),__SIZE_MAX__);
    //thrust::sequence(snps.begin(),snps.begin()+start_size,0);

    //Set parent indexes to be:
    //{junk_slot, 1, 2, 3, 4, 5, ..., start_size+1, __SIZE_MAX__, ...}
    //thrust::fill(parent_indexes.begin()+start_size+1,parent_indexes.end(),__SIZE_MAX__);
    //thrust::sequence(parent_indexes.begin()+1,parent_indexes.begin()+start_size+1,1);
    //thrust::fill(branch_indexes.begin(),branch_indexes.end(),__SIZE_MAX__);

    //Finish parent0 by bitwise NOTing it
    thrust::for_each_n(thrust::make_counting_iterator<std::size_t>(1), start_size, [&](std::size_t index){
        for(std::size_t i = 0; i < snpDataSizeULL-1; ++i)
        {
            parent0[index*snpDataSizeULL+i] =  ~parent0[index*snpDataSizeULL+i];
        }
        parent0[(index+1)*snpDataSizeULL-1] = (~parent0[(index+1)*snpDataSizeULL-1]) & mask;
    });

    //Set offsets to be:
    //{junk_slot, 1, 1, 1, ..., 1, 0, 0, 0, ...}
    //thrust::fill(offsets.begin()+1, offsets.begin()+start_size+1,1);
    //thrust::fill(offsets.begin()+start_size+1, offsets.end(),0);

    //Set direction
    //thrust::fill(direction.begin(),direction.end(),true);

    //Calculate branch0
    thrust::for_each_n(parent_metadata.begin(),
                       num_leaves,
                       [&](ParentMetadata& p){
        std::size_t index = p.snp;
        std::ptrdiff_t offset = p.offset;
        std::size_t parent_index = p.parent_index;
        bool dir = p.direction;
        int count0 = 0, count1 = 0;
        for(std::size_t i = 0; i < snpDataSizeULL-1; ++i)
        {
            unsigned long long val = parent0[parent_index*snpDataSizeULL+i] & ~map_data[(index+offset)*snpDataSizeULL+i];
            count0+=popcount1(val);
            child0[(2*parent_index-2)*snpDataSizeULL+i] = val;
        }
        child0[(2*parent_index-1)*snpDataSizeULL-1] = parent0[(parent_index+1)*snpDataSizeULL-1] & ~map_data[(index+offset+1)*snpDataSizeULL-1] & mask;
        count0+=popcount1(child0[(2*parent_index-1)*snpDataSizeULL-1]);
        for(std::size_t i = 0; i < snpDataSizeULL; ++i)
        {
            unsigned long long val = parent0[parent_index*snpDataSizeULL+i] & map_data[(index+offset)*snpDataSizeULL+i];
            count1+=popcount1(val);
            child0[(2*parent_index-1)*snpDataSizeULL+i] = val;
        }
        //printf("%zd %d %d %d\n", index, offset, count0, count1);
        ChildMetadata c0, c1;
        c0.snp = p.snp;
        c1.snp = p.snp;
        c0.offset = offset+(dir)? 1 : -1;
        c1.offset = offset+(dir)? 1 : -1;
        c0.child_index = 2*parent_index-2;
        c1.child_index = 2*parent_index-1;
        //n*(n-1.0)*0.5
        c0.binom = count0*(count0-1)/2;
        c1.binom = count1*(count1-1)/2;
        bool end = ((offset == 0 && !dir) || (offset == num_snps && dir)) ? true : false;
        c0.hit_end = end;
        c1.hit_end = end;
        nonzero_branches[2*parent_index-2] = (count0) ? 1 : 0;
        nonzero_branches[2*parent_index-1] = (count1) ? 1 : 0;
        child_metadata  [2*parent_index-2] = c0;
        child_metadata  [2*parent_index-1] = c1;
    });

    //Perform scan to get target indexes
    thrust::inclusive_scan(nonzero_branches.begin(), nonzero_branches.begin()+2*num_leaves, compact_branches.begin());

    std::size_t new_leaves = compact_branches[2*num_leaves-1]-1;

    //Compact arrays
    thrust::for_each_n(thrust::make_zip_iterator(thrust::make_tuple(
                           child_metadata.begin(),
                           compact_branches.begin(),
                           nonzero_branches.begin())),
                       2*num_leaves,
                       [&](thrust::tuple<ChildMetadata, std::size_t, std::size_t>&& t){
        ChildMetadata& cm = t.get<0>();
        std::size_t compact_index = t.get<1>();
        std::size_t nonzero = t.get<2>();
        if (nonzero)
        {
            for(std::size_t i = 0; i < snpDataSizeULL; ++i)
            {
                parent0[compact_index*snpDataSizeULL+i] = child0[cm.child_index*snpDataSizeULL+i];
            }
            ParentMetadata pm;
            pm.snp    = cm.snp;
            pm.offset = cm.offset;
            pm.parent_index = compact_index;
            parent_metadata[compact_index-1] = pm;
            binoms[compact_index-1] = cm.binom;
            snps[compact_index-1] = cm.snp;
        }
    });

    //We should now have snps and binoms in a key-value relationship. Now, we add up the binomial coefficients.
    auto new_ends = thrust::reduce_by_key(snps.begin(),
                                          snps.begin()+new_leaves,
                                          binoms.begin(),
                                          reduced_snps.begin(),
                                          reduced_binoms.begin(),
                                          thrust::equal_to<std::size_t>(),
                                          thrust::plus<unsigned long long>());
    std::size_t active_snps = thrust::distance(snps.begin(),new_ends.first);

    //Calculate EHH
    thrust::for_each_n(thrust::make_counting_iterator<std::size_t>(0),
                       active_snps,
                       [&](std::size_t index){
        std::size_t snp = reduced_snps[index];
        ehh0[index] = reduced_binoms[index]/(double)binom_core[snp];
    });

    thrust::set_intersection_by_key(reduced_snps_prev.begin(),reduced_snps_prev.begin()+active_snps_prev, //A keys
                                    reduced_snps.begin(),reduced_snps.begin()+active_snps,                //B keys
                                    ehh0_prev.begin(),                                                    //A vals
                                    reduced_snps_prev_intersection.begin(),                               //Result keys
                                    ehh0_prev_intersection.begin()                                        //Result values
                                    );

    //thrust::set_union_by_key(reduced_snps.begin(),reduced_snps.begin()+active_snps,
    //                              )

    //Compact arrays
    //thrust::for_each_n(thrust::make_zip_iterator(thrust::make_tuple(
    //                       snp_branches.begin()+1,
    //                       offset_branches.begin()+1,
    //                       nonzero_branches.begin()+1,
    //                       branch_counts.begin()+1,
    //                       branch_indexes.begin()+1,
    //                       compact_branches.begin()+1)),
    //                   2*num_leaves,
    //                   [&](thrust::tuple<std::size_t, std::size_t, std::size_t, unsigned int, std::size_t, std::size_t> t){
    //    std::size_t branch_snp = t.get<0>();
    //    std::size_t offset_branches = t.get<1>();
    //    std::size_t nonzero = t.get<2>();
    //    unsigned int count = t.get<3>();
    //    std::size_t branch_index = t.get<4>();
    //    std::size_t compact = t.get<5>();
    //    std::size_t target_index = (nonzero) ? compact : 0;
    //    snps          [target_index] = branch_snp;
    //    offsets       [target_index] = offset_branches;
    //    parent_indexes[target_index] = target_index;
    //    //printf("%zd %zd %zd %zd %d\n", branch_index, nonzero, compact, target_index, count);
    //});

    //calcIhs(hap.value(), map.value(), outfile.value(), cutoff.value(), minMAF.value(), (double) scale.value(), maxExtend.value(), binfac.value(), binom.value());


#if MPI_FOUND
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
#endif
    return ret;
}


