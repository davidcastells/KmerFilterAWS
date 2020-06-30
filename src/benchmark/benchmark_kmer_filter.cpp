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

#include "benchmark_kmer_filter.h"

#include "../filter/kmer_filter.h"
#include "benchmark_edit_alg.h"

/*
 * Benchmark kmer-filter
 * 
 * @param filter_input input parameters
 * @param kmer_length length of the n-grams
 * 
 * for each input we compute the min error bound
 *  The pattern is the short sequence we will be looking into the (bigger) text 
 */
void benchmark_kmer_filter(filter_input_t* const filter_input, const int kmer_length) 
{
  // Parameters
  kmer_counting_nway_t* const kmer_counting = kmer_counting_new(kmer_length, filter_input->mm_allocator);
  
  // Allocate
  // computes the histogram of the pattern
  kmer_counting_pattern_compute_histogram(kmer_counting, (uint8_t*)filter_input->pattern,filter_input->pattern_length);
  
  if (filter_input->verbose) 
  {
      printf("=>Histogram Count\n");
      
      for (int i=0; i < kmer_counting->num_kmers; i++)
      {
          printf("kmer[%d]=%d\n", i, kmer_counting->kmer_count_pattern[i]);
      }  
  }
  
  // Filter
  timer_start(&filter_input->timer);
  const uint64_t min_error_bound = kmer_counting_min_bound(kmer_counting,(uint8_t*)filter_input->text,
          filter_input->text_length,filter_input->max_error);
  
  timer_stop(&filter_input->timer);
  // Check result
  if (filter_input->check) 
  {
    const bool accepted = (min_error_bound <= filter_input->max_error);
    benchmark_check(filter_input,accepted);
  }
  // Free
  kmer_counting_destroy(kmer_counting);
}




