/*
 *  Wavefront Alignments Algorithms
 *  Copyright (c) 2019 by Santiago Marco-Sola  <santiagomsola@gmail.com>
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

#ifndef BENCHMARK_UTILS_H_
#define BENCHMARK_UTILS_H_

#include "../utils/commons.h"
#include "../system/profiler_timer.h"
#include "../system/mm_allocator.h"

/*
 * Filter Input
 */
typedef struct {
  // Sequence ID
  int sequence_id;
  // Pattern
  char* pattern;
  int pattern_length;
  // Text
  char* text;
  int text_length;
  // Error
  int max_error;
  // Profile
  profiler_timer_t timer;
  int candidates_total;
  int candidates_tp;
  int candidates_fp;
  int candidates_tn;
  int candidates_fn;
  // MM
  mm_allocator_t* mm_allocator;
  // DEBUG
  bool check;
  bool verbose;
} filter_input_t;

/*
 * Setup
 */
void filter_input_clear(
    filter_input_t* const filter_input);

#endif /* BENCHMARK_UTILS_H_ */
