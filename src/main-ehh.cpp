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

#include <iostream>

#include "hapbin.hpp"
#include "hapmap.hpp"
#include "ehh.hpp"
#include "argparse.hpp"
#include "ehhfinder.hpp"

int main(int argc, char** argv)
{
    Argument<bool> help('h', "help", "Show this help", true, false);
    Argument<const char*> hap('d', "hap", "Hap file", false, true, "");
    Argument<const char*> map('m', "map", "Map file", false, true, "");
    Argument<double> cutoff('c', "cutoff", "EHH cutoff value (default: 0.05)", false, false, 0.05);
    Argument<double> minMAF('b', "minmaf", "Minimum allele frequency (default: 0.05)", false, false, 0.05);
    Argument<unsigned long long> brTerm('t', "minbranch", "Minimum branch population (default: 1)", false, false, 1ULL);
    Argument<unsigned long long> scale('s', "scale", "Gap scale parameter in bp, used to scale gaps > scale parameter as in Voight, et al.", false, false, 20000ULL);
    Argument<const char*> locus('l', "locus", "Locus", false, true, 0);
    ArgParse argparse({&help, &hap, &map, &locus, &cutoff, &minMAF, &scale, &brTerm}, "Usage: ehhbin --map input.map --hap input.hap --locus id");
    argparse.parseArguments(argc, argv);
    //using HapMapType = HapMap<CTCBitset<2*INDIVIDUALS>>;
    using HapMapType = HapMap;
    HapMapType hmap;
    if (!hmap.loadHap(hap.value()))
    {
        return 1;
    }
    hmap.loadMap(map.value());
    EHH e;
    std::size_t l = hmap.idToLine(locus.value());
    if (l == std::numeric_limits<std::size_t>::max())
    {
        std::cerr << "no locus with the id: " << locus.value() << std::endl;
        return 2;
    }
    EHHFinder finder(hmap.snpDataSize(), hmap.snpDataSize(), 1000, cutoff.value(), minMAF.value(), (double) scale.value(), brTerm.value());
    e = finder.find(&hmap, l, true);
    e.printEHH(&hmap);
    std::cout << "iHS: " << log(e.iHH_a/e.iHH_d) << std::endl;
    std::cout << "MAF: " << (double)e.num/(double)hmap.snpLength() << std::endl;
}
