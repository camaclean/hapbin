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

#include "ehhfinder.hpp"
#include "hapmap.hpp"
#include <algorithm>

EHHFinder::EHHFinder(HapMap* hmA, HapMap* hmB, std::size_t maxBreadth, double cutoff, double minMAF, double scale, double brCutoff)
    : m_hmA(hmA)
    , m_hmB(hmB)
    , m_maxBreadth(maxBreadth)
    , m_cutoff(cutoff)
    , m_minMAF(minMAF)
    , m_scale(scale)
    , m_freqA{}
    , m_freqB{}
    , m_freqP{}
    , m_ehhA{}
    , m_ehhB{}
    , m_ehhP{}
    , m_brCutoff{brCutoff}
    , m_term0counts(NULL)
    , m_term1counts(NULL)
    , m_termPooledCounts(NULL)
    , m_term0countsTmp(NULL)
    , m_term1countsTmp(NULL)
    , m_termPooledCountsTmp(NULL)
{
    if (m_hmB == NULL) {
        m_hdA = hmA->rawData();
        m_snpDataSizeA = hmA->snpDataSize();
        m_snpDataSizeULL_A = hmA->snpDataSizeULL();
        m_snpLengthA = hmA->snpLength();
#if VEC==4
        m_maskA = ::bitsetMask4(hmA->snpLength());
#elif VEC==2
        m_maskA = ::bitsetMask2(hmA->snpLength());
#else
        m_maskA = ::bitsetMask<HapMap::PrimitiveType>(hmA->snpLength());
#endif
        m_brTermMax = m_brCutoff * m_hmA->snpLength();
        m_parent0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeA*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
        m_parent1 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeA*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
        m_branch0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeA*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
        m_branch1 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeA*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
    }
    else
    {
        m_hdA = hmA->rawData();
        m_snpDataSizeA = hmA->snpDataSize();
        m_snpDataSizeULL_A = hmA->snpDataSizeULL();
        m_hdB = hmB->rawData();
        m_snpDataSizeB = hmB->snpDataSize();
        m_snpDataSizeULL_B = hmB->snpDataSizeULL();
        m_snpLengthA = hmA->snpLength();
        m_snpLengthB = hmB->snpLength();
#if VEC==4
        m_maskA = ::bitsetMask4(hmA->snpLength());
        m_maskB = ::bitsetMask4(hmB->snpLength());
#elif VEC==2
        m_maskA = ::bitsetMask2(hmA->snpLength());
        m_maskB = ::bitsetMask2(hmB->snpLength());
#else
        m_maskA = ::bitsetMask<HapMap::PrimitiveType>(hmA->snpLength());
        m_maskB = ::bitsetMask<HapMap::PrimitiveType>(hmB->snpLength());
#endif
        m_brTermMax = m_brCutoff * (m_hmA->snpLength() + m_hmB->snpLength());
        m_parent0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeA*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
        m_parent1 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeB*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
        m_branch0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeA*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
        m_branch1 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeB*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
    }
    if (m_brTermMax == 0) m_brTermMax = 1;
    m_term0counts = new size_t[m_brTermMax]();
    m_term1counts = new size_t[m_brTermMax]();
    m_termPooledCounts = new size_t[m_brTermMax]();
    m_term0countsTmp = new size_t[m_brTermMax]();
    m_term1countsTmp = new size_t[m_brTermMax]();
    m_termPooledCountsTmp = new size_t[m_brTermMax]();
    
}

/**
 * Set initial state. Set m_parent0 to '0' core haplotype positions, m_parent1 to '1' core haplotype positions.
 * 
 * calcBranch counts the state of the parent (previous) branch, not the counts of the new (current) level. Therefore, we must advance
 * the calculations by one iteration before starting. 
 */
void EHHFinder::setInitial(std::size_t focus, std::size_t line)
{
    m_parent0count = 2ULL;
    m_parent1count = 2ULL;
    m_branch0count = 0ULL;
    m_branch1count = 0ULL;
    std::fill_n(m_term0counts, m_brTerm0, 0ULL);
    std::fill_n(m_term1counts, m_brTerm1, 0ULL);

    for (int j = 0; j < m_snpDataSizeA; ++j)
    {
        m_parent0[               j] =  ~m_hdA[focus*m_snpDataSizeA+j]  &   m_hdA[line*m_snpDataSizeA+j];
        m_parent0[m_snpDataSizeA+j] = (~m_hdA[focus*m_snpDataSizeA+j]) & (~m_hdA[line*m_snpDataSizeA+j]);
        m_parent1[               j] =   m_hdA[focus*m_snpDataSizeA+j]  &   m_hdA[line*m_snpDataSizeA+j];
        m_parent1[m_snpDataSizeA+j] =   m_hdA[focus*m_snpDataSizeA+j]  &  ~m_hdA[line*m_snpDataSizeA+j];
    }
    m_parent0[  m_snpDataSizeA-1] &= m_maskA;
    m_parent0[2*m_snpDataSizeA-1] &= m_maskA;
    m_parent1[2*m_snpDataSizeA-1] &= m_maskA;
}

void EHHFinder::setInitialXPEHH(std::size_t focus)
{
    for (std::size_t j = 0; j < m_snpDataSizeA; ++j)
    {
        m_parent0[j] = ~m_hdA[focus*m_snpDataSizeA+j];
    }
    m_parent0[m_snpDataSizeA-1] &= m_maskA;
    for (std::size_t j = 0; j < m_snpDataSizeA; ++j)
    {
        m_parent0[j+m_snpDataSizeA] = m_hdA[focus*m_snpDataSizeA+j];
    }
    for (std::size_t j = 0; j < m_snpDataSizeB; ++j)
    {
        m_parent1[j] = ~m_hdB[focus*m_snpDataSizeB+j];
    }
    m_parent1[m_snpDataSizeB-1] &= m_maskB;
    for (std::size_t j = 0; j < m_snpDataSizeB; ++j)
    {
        m_parent1[j+m_snpDataSizeB] = m_hdB[focus*m_snpDataSizeB+j];
    }
}

void EHHFinder::calcBranch(std::size_t focus, HapMap::PrimitiveType* parent, std::size_t parentcount, HapMap::PrimitiveType* branch, std::size_t& branchcount, std::size_t currLine, double freq, double &probs, std::size_t* termCts, std::size_t brTerm, bool* overflow)
{
    std::size_t bcnt = 0;
    for(std::size_t i = 0; i < parentcount; ++i)
    {
        int count = 0;
        unsigned long long *leaf = (unsigned long long*) &parent[i*m_snpDataSizeA];
        for (std::size_t j = 0; j < m_snpDataSizeULL_A; ++j)
        {
            count += popcount1(leaf[j]);
        }
        
        if (count == 0)
        {
            continue;
        }
        else if (count <= brTerm)
        {
            ++termCts[count-1];
            continue;
        }
        else
        {
            probs += (count*freq)*(count*freq);
            for(std::size_t j = 0; j < m_snpDataSizeA; ++j)
            {
                branch[bcnt*m_snpDataSizeA+j] = parent[i*m_snpDataSizeA+j] & m_hdA[currLine*m_snpDataSizeA+j];
            }
            ++bcnt;
            for(std::size_t j = 0; j < m_snpDataSizeA-1; ++j)
            {
                branch[bcnt*m_snpDataSizeA+j] = parent[i*m_snpDataSizeA+j] & ~m_hdA[currLine*m_snpDataSizeA+j];
            }
            branch[bcnt*m_snpDataSizeA+m_snpDataSizeA-1] = (parent[i*m_snpDataSizeA+m_snpDataSizeA-1] & ~m_hdA[currLine*m_snpDataSizeA+m_snpDataSizeA-1]) & m_maskA;
            ++bcnt;
        }
        if (bcnt > m_maxBreadth-2)
        {
            *overflow = true;
            return;
        }
    }
    branchcount = bcnt;
}

void EHHFinder::calcBranchXPEHH(std::size_t currLine, bool* overflow)
{
    for(std::size_t i = 0; i < m_parent0count; ++i)
    {
        int numA = 0;
        unsigned long long *parentA_ULL = (unsigned long long*) &m_parent0[i*m_snpDataSizeA];
        for (std::size_t j = 0; j < m_snpDataSizeULL_A; ++j)
        {
            numA += popcount1(parentA_ULL[j]);
        }
        int numB = 0;
        unsigned long long *parentB_ULL = (unsigned long long*) &m_parent1[i*m_snpDataSizeB];
        for (size_t j = 0; j < m_snpDataSizeULL_B; ++j)
        {
            numB += popcount1(parentB_ULL[j]);
        }
        int numPooled = numA + numB;
        //std::cout << i << "numA: " << numA << " numB: " << numB << " numP: " << numPooled << " freqP: " << m_freqP << " " << m_snpDataSizeUL_A  << std::endl;
        if (numPooled == 0)
        {
            continue;
        } 
        else if (numPooled <= ( m_brTermMax))
        {
            if (numA != 0)
                ++m_term0countsTmp[numA-1];
            if (numB != 0)
                ++m_term1countsTmp[numB-1];
            continue;
        } 
        else
        {
            m_ehhA += (numA*m_freqA)*(numA*m_freqA);
            m_ehhB += (numB*m_freqB)*(numB*m_freqB);
            m_ehhP += (numPooled*m_freqP)*(numPooled*m_freqP);
            
            //Population A
            for(std::size_t j = 0; j < m_snpDataSizeA; ++j)
            {
                m_branch0[m_branch0count*m_snpDataSizeA+j] = m_parent0[i*m_snpDataSizeA+j] & m_hdA[currLine*m_snpDataSizeA+j];
            }
            ++m_branch0count;
            for(std::size_t j = 0; j < m_snpDataSizeA-1; ++j)
            {
                m_branch0[m_branch0count*m_snpDataSizeA+j] = m_parent0[i*m_snpDataSizeA+j] & ~m_hdA[currLine*m_snpDataSizeA+j];
            }
            m_branch0[m_branch0count*m_snpDataSizeA+m_snpDataSizeA-1] = (m_parent0[i*m_snpDataSizeA+m_snpDataSizeA-1] & ~m_hdA[currLine*m_snpDataSizeA+m_snpDataSizeA-1]) & m_maskA;
            ++m_branch0count;
            
            //Population B
            for(std::size_t j = 0; j < m_snpDataSizeB; ++j)
            {
                m_branch1[m_branch1count*m_snpDataSizeB+j] = m_parent1[i*m_snpDataSizeB+j] & m_hdB[currLine*m_snpDataSizeB+j];
            }
            ++m_branch1count;
            for(std::size_t j = 0; j < m_snpDataSizeB-1; ++j)
            {
                m_branch1[m_branch1count*m_snpDataSizeB+j] = m_parent1[i*m_snpDataSizeB+j] & ~m_hdB[currLine*m_snpDataSizeB+j];
            }
            m_branch1[m_branch1count*m_snpDataSizeB+m_snpDataSizeB-1] = (m_parent1[i*m_snpDataSizeB+m_snpDataSizeB-1] & ~m_hdB[currLine*m_snpDataSizeB+m_snpDataSizeB-1]) & m_maskB;
            ++m_branch1count;
        }
        if (m_branch0count > m_maxBreadth-2)
        {
            *overflow = true;
            return;
        }
    }
}

void EHHFinder::calcBranchesXPEHH(std::size_t currLine)
{
    bool overflow;
    bool realloced = false;
    m_ehhA = 0.0;
    m_ehhB = 0.0;
    m_ehhP = 0.0;
    while (true)
    {
        overflow = false;
        std::fill_n(m_term0countsTmp,m_brTermMax, 0ULL);
        std::fill_n(m_term1countsTmp,m_brTermMax, 0ULL);
        m_branch0count = 0;
        m_branch1count = 0;
        calcBranchXPEHH(currLine, &overflow);
        if (overflow)
        {
            m_maxBreadth += 500;
            aligned_free(m_branch0);
            aligned_free(m_branch1);
            m_branch0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeA*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
            m_branch1 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeB*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
            realloced = true;
        }
        else
            break;
    }
    if (realloced)
    {
        aligned_free(m_parent0);
        aligned_free(m_parent1);
        m_parent0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeA*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
        m_parent1 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeB*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
    }
    for(std::size_t i = 0; i < m_brTermMax; ++i)
    {
        m_term0counts[i] += m_term0countsTmp[i];
    }
    for(std::size_t i = 0; i < m_brTermMax; ++i)
    {
        m_term1counts[i] += m_term1countsTmp[i];
    }
    m_parent0count = m_branch0count;
    m_parent1count = m_branch1count;
    m_branch0count = 0ULL;
    m_branch1count = 0ULL;
    std::swap(m_parent0, m_branch0);
    std::swap(m_parent1, m_branch1);
}

std::pair< EHH, EHH > EHHFinder::findXPEHH(std::size_t focus, std::atomic<unsigned long long>* reachedEnd)
{
    m_parent0count = 2ULL;
    m_parent1count = 2ULL;
    m_branch0count = 0ULL;
    m_branch1count = 0ULL;
    std::fill_n(m_term0counts,m_brTermMax, 0ULL);
    std::fill_n(m_term1counts,m_brTermMax, 0ULL);
    std::pair<EHH,EHH> ret;
    ret.first.index = focus;
    for(std::size_t i = 0; i < m_snpDataSizeA; ++i)
    {
        ret.first.num += POPCOUNT(m_hdA[focus*m_snpDataSizeA+i]);
    }
    ret.first.numNot = m_snpLengthA - ret.first.num;
    for(std::size_t i = 0; i < m_snpDataSizeB; ++i)
    {
        ret.second.num += POPCOUNT(m_hdB[focus*m_snpDataSizeB+i]);
    }
    ret.second.numNot = m_snpLengthB - ret.second.num;
    m_freqA = 1.0/(double)m_snpLengthA;
    m_freqB = 1.0/(double)m_snpLengthB;
    m_freqP = 1.0/(double)(m_snpLengthA + m_snpLengthB);
    double probASingle = m_freqA*m_freqA;
    double probBSingle = m_freqB*m_freqB;
    double probPSingle = m_freqP*m_freqP;
    double f = ret.first.num*m_freqA;
    double lastEhhA = f*f+(1.0-f)*(1.0-f);
    f = ret.second.num*m_freqB;
    double lastEhhB = f*f+(1.0-f)*(1.0-f);
    f = (ret.first.num+ret.second.num)*m_freqP;
    double lastEhhP = f*f+(1.0-f)*(1.0-f);
    
    setInitialXPEHH(focus);
    calcBranchesXPEHH(focus-1);
    std::fill_n(m_term0counts,m_brTermMax, 0ULL);
    std::fill_n(m_term1counts,m_brTermMax, 0ULL);
    if (focus > 1)
    {
        for (std::size_t currLine = focus - 2;; --currLine)
        {
            
            double scale = (double)(m_scale) / (double)(m_hmA->physicalPosition(currLine+2) - m_hmA->physicalPosition(currLine+1));
            if (scale > 1)
                scale=1;
            
            calcBranchesXPEHH(currLine);
            
            for(std::size_t i = 1; i <= m_brTermMax; ++i)
            {
                m_ehhA += i*i*m_term0counts[i-1]*probASingle;
                m_ehhB += i*i*m_term1counts[i-1]*probBSingle;
                m_ehhP += i*i*(m_term0counts[i-1]+m_term1counts[i-1])*probPSingle;
            }
            
            if (lastEhhP <= m_cutoff + 1e-15)
                break;
            ret.first.iHH_d  += (m_hmA->geneticPosition(currLine+2)-m_hmA->geneticPosition(currLine+1))*(lastEhhA + m_ehhA)*scale*0.5;    
            ret.second.iHH_d += (m_hmA->geneticPosition(currLine+2)-m_hmA->geneticPosition(currLine+1))*(lastEhhB + m_ehhB)*scale*0.5;
            
            lastEhhA = m_ehhA;
            lastEhhB = m_ehhB;
            lastEhhP = m_ehhP;
            
            if (m_parent0count == 0 && m_parent1count == 0)
                break;
            if(currLine == 0)
            {
                ++(*reachedEnd);
                return std::pair<EHH,EHH>();
            }
        }
    }
    
    f = ret.first.num*m_freqA;
    lastEhhA = f*f+(1.0-f)*(1.0-f);
    f = ret.second.num*m_freqB;
    lastEhhB = f*f+(1.0-f)*(1.0-f);
    f = (ret.first.num+ret.second.num)*m_freqP;
    lastEhhP = f*f+(1.0-f)*(1.0-f);
    
    setInitialXPEHH(focus);
    calcBranchesXPEHH(focus+1);
    std::fill_n(m_term0counts, m_brTermMax, 0ULL);
    std::fill_n(m_term1counts, m_brTermMax, 0ULL);
    for (std::size_t currLine = focus + 2; currLine < m_hmA->numSnps(); ++currLine)
    {
        double scale = (double)(m_scale) / (double)(m_hmA->physicalPosition(currLine-1) - m_hmA->physicalPosition(currLine-2));
        if (scale > 1)
            scale=1;
        
        calcBranchesXPEHH(currLine);
        
        for(std::size_t i = 1; i <= m_brTermMax; ++i)
        {
            m_ehhA += i*i*m_term0counts[i-1]*probASingle;
            m_ehhB += i*i*m_term1counts[i-1]*probBSingle;
            m_ehhP += i*i*(m_term0counts[i-1]+m_term1counts[i-1])*probPSingle;
        }
        
        if (lastEhhP <= m_cutoff + 1e-15)
            break;
        ret.first.iHH_d  += (m_hmA->geneticPosition(currLine-1)-m_hmA->geneticPosition(currLine-2))*(lastEhhA + m_ehhA)*scale*0.5;    
        ret.second.iHH_d += (m_hmA->geneticPosition(currLine-1)-m_hmA->geneticPosition(currLine-2))*(lastEhhB + m_ehhB)*scale*0.5;
        
        lastEhhA = m_ehhA;
        lastEhhB = m_ehhB;
        lastEhhP = m_ehhP;
        
        if (m_parent0count == 0 && m_parent1count == 0)
            break;
        if (currLine == m_hmA->numSnps()-1)
        {
            ++(*reachedEnd);
            return std::pair<EHH,EHH>();
        }
    }
    return ret;
}

void EHHFinder::calcBranches(std::size_t focus, std::size_t currLine, double freq0,  double freq1, HapStats &stats)
{
    bool overflow;
    bool realloced = false;
    while(true)
    {
        overflow = false;
        std::fill_n(m_term0countsTmp, m_brTerm0, 0ULL);
        std::fill_n(m_term1countsTmp, m_brTerm1, 0ULL);
        calcBranch(focus, m_parent0, m_parent0count, m_branch0, m_branch0count, currLine, freq0, stats.probsNot, m_term0countsTmp, m_brTerm0, &overflow);
        if(overflow)
        {
            m_maxBreadth += 100;
            aligned_free(m_branch0);
            aligned_free(m_branch1);
            m_branch0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeA*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
            m_branch1 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeB*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
            realloced = true;
            continue;
        }
        overflow = false;
        calcBranch(focus, m_parent1, m_parent1count, m_branch1, m_branch1count, currLine, freq1, stats.probs, m_term1countsTmp, m_brTerm1, &overflow);
        if(overflow)
        {
            m_maxBreadth += 100;
            aligned_free(m_branch0);
            aligned_free(m_branch1);
            m_branch0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeA*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
            m_branch1 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeB*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
            realloced = true;
            continue;
        }
        break;
    }
    if (realloced)
    {
        aligned_free(m_parent0);
        aligned_free(m_parent1);
        m_parent0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeA*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
        m_parent1 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeB*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
    }
    m_parent0count = m_branch0count;
    m_parent1count = m_branch1count;
    m_branch0count = 0ULL;
    m_branch1count = 0ULL;
    for(std::size_t i = 0; i < m_brTerm0; ++i)
    {
        m_term0counts[i] += m_term0countsTmp[i];
    }
    for(std::size_t i = 0; i < m_brTerm1; ++i)
    {
        m_term1counts[i] += m_term1countsTmp[i];
    }
    std::swap(m_parent0, m_branch0);
    std::swap(m_parent1, m_branch1);
}

EHH EHHFinder::find(std::size_t focus, std::atomic<unsigned long long>* reachedEnd, std::atomic<unsigned long long>* outsideMaf, bool ehhsave)
{
    m_parent0count = 2ULL;
    m_parent1count = 2ULL;
    m_branch0count = 0ULL;
    m_branch1count = 0ULL;
    EHH ret;
    ret.index = focus;
    
    for(std::size_t i = 0; i < m_snpDataSizeA; ++i)
    {
        ret.num += POPCOUNT(m_hdA[focus*m_snpDataSizeA+i]);
    }
    ret.numNot = m_snpLengthA - ret.num;
    
    m_brTerm0 = ret.numNot*m_brCutoff;
    m_brTerm1 = ret.num*m_brCutoff;
    if (m_brTerm0 == 0)
        m_brTerm0 = 1;
    if (m_brTerm1 == 0)
        m_brTerm1 = 1;

    double maxEHH = ret.num/(double)m_snpLengthA;
    if (!(maxEHH <= 1.0 - m_minMAF && maxEHH >= m_minMAF))
    {
        ++(*outsideMaf);
        return EHH();
    }
    if (focus < 2)
    {
        ++(*reachedEnd);
        return EHH();
    }
    
    double freq0 = 1.0/(double)ret.numNot;
    double freq1 = 1.0/(double)ret.num;
    double probSingle = freq1*freq1;
    double probNotSingle = freq0*freq0;
    double probs, probsNot;
    double lastProbs = 1.0, lastProbsNot = 1.0;
    
    setInitial(focus, focus-1);
    
    for (std::size_t currLine = focus - 2;; --currLine)
    {
        HapStats stats;
        double scale = (double)(m_scale) / (double)(m_hmA->physicalPosition(currLine+2) - m_hmA->physicalPosition(currLine+1));
        if (scale > 1)
            scale=1;
        
        calcBranches(focus, currLine, freq0, freq1, stats);
        
        for(std::size_t i = 1; i <= m_brTerm1; ++i)
        {
            stats.probs += i*i*m_term1counts[i-1]*probSingle;
        }
        for(std::size_t i = 1; i <= m_brTerm0; ++i)
        {
            stats.probsNot += i*i*m_term0counts[i-1]*probNotSingle;
        }

        if (lastProbs > m_cutoff + 1e-15)
            ret.iHH_d += (m_hmA->geneticPosition(currLine+2)-m_hmA->geneticPosition(currLine+1))*(lastProbs + stats.probs)*scale*0.5;    
        if (lastProbsNot > m_cutoff + 1e-15)
            ret.iHH_a += (m_hmA->geneticPosition(currLine+2)-m_hmA->geneticPosition(currLine+1))*(lastProbsNot + stats.probsNot)*scale*0.5;   
        
        lastProbs = stats.probs;
        lastProbsNot = stats.probsNot;
        if (ehhsave)
            ret.upstream.push_back(std::move(stats));
        
        if (lastProbs <= m_cutoff + 1e-15 && lastProbsNot <= m_cutoff + 1e-15)
            break;
        if (m_parent0count == 0 && m_parent1count == 0)
            break;
        if (currLine == 0)
        {
            ++(*reachedEnd);
            return EHH();
        }
    }
    
    setInitial(focus,focus+1);
    lastProbs = 1.0, lastProbsNot = 1.0;
    for (std::size_t currLine = focus + 2; currLine < m_hmA->numSnps(); ++currLine)
    {
        HapStats stats;
        double scale = double(m_scale) / double(m_hmA->physicalPosition(currLine-1) - m_hmA->physicalPosition(currLine-2));
        if (scale > 1)
            scale=1;
        
        int core0 = 0, core1 = 0;
        
        calcBranches(focus, currLine, freq0, freq1, stats);
        
        for(std::size_t i = 1; i <= m_brTerm1; ++i)
        {
            stats.probs += i*i*m_term1counts[i-1]*probSingle;
        }
        for(std::size_t i = 1; i <= m_brTerm0; ++i)
        {
            stats.probsNot += i*i*m_term0counts[i-1]*probNotSingle;
        }
        
        if (lastProbs > m_cutoff + 1e-15) {
            ret.iHH_d += (m_hmA->geneticPosition(currLine-1)-m_hmA->geneticPosition(currLine-2))*(lastProbs + stats.probs)*scale*0.5;
        }
        if (lastProbsNot > m_cutoff + 1e-15) {
            ret.iHH_a += (m_hmA->geneticPosition(currLine-1)-m_hmA->geneticPosition(currLine-2))*(lastProbsNot + stats.probsNot)*scale*0.5;
        }
        
        lastProbs = stats.probs;
        lastProbsNot = stats.probsNot;
        if (ehhsave)
            ret.downstream.push_back(std::move(stats));
        
        if (lastProbs <= m_cutoff + 1e-15 && lastProbsNot <= m_cutoff + 1e-15)
            break;
        if (m_parent0count == 0 && m_parent1count == 0)
            break;
        
        if (currLine == m_hmA->numSnps()-1)
        {
            ++(*reachedEnd);
            return EHH();
        }
    }
    return ret;
}

EHHFinder::~EHHFinder()
{
    aligned_free(m_branch0);
    aligned_free(m_branch1);
    aligned_free(m_parent0);
    aligned_free(m_parent1);
    delete m_term0counts;
    delete m_term1counts;
    delete m_termPooledCounts;
    delete m_term0countsTmp;
    delete m_term1countsTmp;
    delete m_termPooledCountsTmp;
}
