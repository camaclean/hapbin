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
#include <limits>

int main(int argc, char** argv)
{
#if MPI_FOUND
    MPI_Init(&argc, &argv);
#endif
    Argument<bool> help('h', "help", "Show this help", true, false);
    Argument<std::string> hapA('d', "hapA", "Hap file for population A", false, true, "");
    Argument<std::string> hapB('e', "hapB", "Hap file for population B", false, true, "");
    Argument<std::string> map('m', "map", "Map file", false, true, "");
    Argument<double> cutoff('c', "cutoff", "EHH cutoff value (default: 0.05)", false, false, 0.05);
    Argument<double> minMAF('f', "minmaf", "Minimum allele frequency (default: 0.05)", false, false, 0.00);
    Argument<double> binfac('b', "bin", "Frequency bin size (default: 0.02)", false, false, 0.02);
    Argument<unsigned long long> brTerm('t', "minbranch", "Minimum branch population (default: 1)", false, false, 1ULL);
    Argument<unsigned long long> scale('s', "scale", "Gap scale parameter in bp, used to scale gaps > scale parameter as in Voight, et al.", false, false, 20000);
    Argument<std::string> outfile('o', "out", "Output file", false, false, "out.txt");
    ArgParse argparse({&help, &hapA, &hapB, &map, &outfile, &cutoff, &minMAF, &scale, &binfac, &brTerm}, "Usage: xpehhbin --map input.map --hapA inputA.hap --hapB inputB.hap");
    if (!argparse.parseArguments(argc, argv))
        return 1;
    
    std::size_t numSnps;
    numSnps = HapMap::querySnpLength(hapA.value().c_str());
    std::cout << "Haplotypes in population A: " << numSnps << std::endl;
    numSnps = HapMap::querySnpLength(hapB.value().c_str());
    std::cout << "Haplotypes in population B: " << numSnps << std::endl;
    
    calcXpehh(hapA.value(), hapB.value(), map.value(), outfile.value(), cutoff.value(), minMAF.value(), (double) scale.value(), binfac.value(), brTerm.value());
#if MPI_FOUND
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
#endif
    return 0;
}


