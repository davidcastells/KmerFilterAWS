/*
 * Copyright (C) 2020 Universitat Autonoma de Barcelona - David Castells-Rufas <david.castells@uab.cat>
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
 */


/* 
 * File:   FPGAKmerFilter.h
 * Author: dcr
 *
 * Created on February 20, 2020, 6:04 PM
 */

#ifndef FPGAKMERFILTER_H
#define FPGAKMERFILTER_H

#include "OpenCLUtils.h"
#include "../benchmark/benchmark_utils.h"

#include <string>
#include <vector>

using namespace std;

class FPGAKmerFilter 
{
public:
    FPGAKmerFilter();
    virtual ~FPGAKmerFilter();
    
public:
    void initOpenCL(int platform_id);
    void initKernels(int version, string openCLKernelType);
    void finalizeOpenCL();
    void finalizeKernels();
    void addInput(filter_input_t* const filter_input, const int kmer_length);
    void computeAll();
    void destroy();
    void invokeKernel(unsigned char* pattern, unsigned int patternSize, unsigned int* patternIdx,
                                    unsigned char* text, unsigned int textSize, unsigned int* textIdx,
                                    unsigned int* workload, unsigned int tasks);
    
private:
    size_t countRequiredMemory(vector<int>& len);
    void encodeSequence(string bases, unsigned int basesLength, unsigned char* pattern, unsigned int offset);
    
public:
    bool m_verbose;

private:
    vector<string> m_basesPattern;
    vector<string> m_basesText;
    vector<int> m_basesPatternLength;
    vector<int> m_basesTextLength;
    vector<filter_input_t*> m_original;   

    cl_platform_id m_platform;
    cl_context m_context;
    cl_device_id m_deviceId;
    cl_program m_program;
    cl_command_queue m_queue;
    cl_kernel m_kmerKernel;
    string m_openCLFilesPath;
    
    cl_mem m_memPattern;
    cl_mem m_memPatternIdx;
    cl_mem m_memText;
    cl_mem m_memTextIdx;
    cl_mem m_memWorkload;
    
    int m_openCLKernelVersion;
};

#endif /* FPGAKMERFILTER_H */

