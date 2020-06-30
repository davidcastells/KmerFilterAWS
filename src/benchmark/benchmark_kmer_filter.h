/*
 *  Wavefront Alignments Algorithms
 *  Copyright (c) 2020 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of Wavefront Alignments Algorithms.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: Fast Mapping-Candidates Filtering Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#ifndef BENCHMARK_KMER_FILTER_H_
#define BENCHMARK_KMER_FILTER_H_

#include "../utils/commons.h"
#include "../benchmark/benchmark_utils.h"

/*
 * Benchmark kmer-filter
 */
void benchmark_kmer_filter(
    filter_input_t* const filter_input,
    const int kmer_length);

#endif /* BENCHMARK_KMER_FILTER_H_ */
