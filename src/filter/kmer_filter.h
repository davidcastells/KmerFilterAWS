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

#ifndef KMER_FILTER_NWAY_H_
#define KMER_FILTER_NWAY_H_

#include "../utils/commons.h"
#include "../system/mm_allocator.h"

/*
 * Kmer counting filter
 */
typedef uint16_t kmer_count_int_t;        // Counter size
typedef struct {
  // Filter parameters
  uint64_t kmer_length;                   // Kmer length
  uint64_t kmer_mask;                     // Kmer mask to extract kmer offset
  uint64_t num_kmers;                     // Total number of possible kmers in table
  // Key
  uint8_t* key;                           // Key
  uint64_t key_length;                    // Key length
  uint64_t num_key_kmers;                 // Total kmers in chunk
  // Text
  uint64_t curr_text_kmers;               // Current number of kmers contained in text (wrt pattern profile)
  uint64_t max_text_kmers;                // Maximum number of kmers contained in text (wrt pattern profile)
  // Profile tables
  kmer_count_int_t* kmer_count_text;      // Text profile (kmers on text)
  kmer_count_int_t* kmer_count_pattern;   // Key chunks profile (kmers on each key chunk)
  // MM
  mm_allocator_t* mm_allocator;           // MM-Allocator
} kmer_counting_nway_t;

/*
 * Setup
 */
kmer_counting_nway_t* kmer_counting_new(
    const uint64_t kmer_length,
    mm_allocator_t* const mm_allocator);
void kmer_counting_destroy(
    kmer_counting_nway_t* const kmer_counting);

/*
 * Compile Pattern
 */
void kmer_counting_pattern_compute_histogram(
    kmer_counting_nway_t* const kmer_counting,
    uint8_t* const key,
    const uint64_t key_length);

/*
 * Kmer-filter (Compute minimum error bound)
 */
uint64_t kmer_counting_min_bound(
    kmer_counting_nway_t* const kmer_counting,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t max_error);

#endif /* KMER_FILTER_NWAY_H_ */
