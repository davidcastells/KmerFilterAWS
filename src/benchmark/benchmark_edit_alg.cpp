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

#include "benchmark_edit_alg.h"
#include "../alignment/edit_bpm_distance.h"
#include "../alignment/edit_dp.h"

#include <string>

using namespace std;

/*
 * Check
 */
void benchmark_check(
    filter_input_t* const filter_input,
    const bool accepted) {
  // Parameters
  edit_column_t edit_column;
  // Compute Edit Distance
  edit_column_allocate(
      &edit_column,filter_input->pattern_length,
      filter_input->text_length,filter_input->mm_allocator);
  const int edit_distance = edit_dp_distance(&edit_column,
      filter_input->pattern,filter_input->pattern_length,
      filter_input->text,filter_input->text_length);
  // Check result
  ++(filter_input->candidates_total);
  if (accepted) { // It was accepted
    if (edit_distance <= filter_input->max_error) {
      ++(filter_input->candidates_tp);
    } else {
      ++(filter_input->candidates_fp);
    }
  } else { // It was discarded
    if (edit_distance <= filter_input->max_error) {
      ++(filter_input->candidates_fn);
    } else {
      ++(filter_input->candidates_tn);
    }
  }
  // Free
  edit_column_free(&edit_column,filter_input->mm_allocator);
}

/*
 * My Check
 */
void my_benchmark_check(filter_input_t* const filter_input,
	string pattern,
	string text, 
const bool accepted)
{
  // Parameters
  edit_column_t edit_column;
  // Compute Edit Distance
  edit_column_allocate(&edit_column, pattern.size(), text.size(),filter_input->mm_allocator);

  const int edit_distance = edit_dp_distance(&edit_column,
      pattern.c_str(), pattern.size(),
      text.c_str(),text.size());

  // Check result
  ++(filter_input->candidates_total);
  if (accepted) 
  { 
    // It was accepted
    if (edit_distance <= filter_input->max_error) {
      ++(filter_input->candidates_tp);
    } else {
      ++(filter_input->candidates_fp);
    }
  } else { // It was discarded
    if (edit_distance <= filter_input->max_error) {
      ++(filter_input->candidates_fn);
    } else {
      ++(filter_input->candidates_tn);
    }
  }
  // Free
  edit_column_free(&edit_column,filter_input->mm_allocator);
}


void benchmark_check_edit_distance(
    filter_input_t* const filter_input,
    const int edit_distance) {
  benchmark_check(filter_input,edit_distance<=filter_input->max_error);
}
/*
 * Benchmark Edit
 */
void benchmark_edit_dp(
    filter_input_t* const filter_input,
    const int bandwidth) {
  // Parameters
  edit_column_t edit_column;
  // Allocate
  edit_column_allocate(
      &edit_column,filter_input->pattern_length,
      filter_input->text_length,filter_input->mm_allocator);
  // Align
  int edit_distance;
  if (bandwidth == -1) {
    timer_start(&filter_input->timer);
    edit_distance = edit_dp_distance(&edit_column,
        filter_input->pattern,filter_input->pattern_length,
        filter_input->text,filter_input->text_length);
    timer_stop(&filter_input->timer);
  } else {
    timer_start(&filter_input->timer);
    edit_distance = edit_dp_distance_banded(&edit_column,
        filter_input->pattern,filter_input->pattern_length,
        filter_input->text,filter_input->text_length,bandwidth);
    timer_stop(&filter_input->timer);
  }
  // Check result
  ++(filter_input->candidates_total);
  if (edit_distance <= filter_input->max_error) {
    ++(filter_input->candidates_tp);
  } else {
    ++(filter_input->candidates_tn);
  }
  // Free
  edit_column_free(&edit_column,filter_input->mm_allocator);
}
void benchmark_edit_bpm(
    filter_input_t* const filter_input,
    const int bandwidth) {
  // Allocate
  bpm_pattern_t bpm_pattern;
  edit_bpm_pattern_compile(
      &bpm_pattern,filter_input->pattern,
      filter_input->pattern_length,filter_input->mm_allocator);
  // Align
  int edit_distance;
  if (bandwidth == -1) {
    timer_start(&filter_input->timer);
    edit_distance = edit_bpm_distance_compute(
        &bpm_pattern,filter_input->text,filter_input->text_length);
    timer_stop(&filter_input->timer);
  } else {
    timer_start(&filter_input->timer);
    edit_distance = edit_bpm_distance_compute_cutoff(
          &bpm_pattern,filter_input->text,
          filter_input->text_length,bandwidth,true);
    timer_stop(&filter_input->timer);
  }
  // Check result
  ++(filter_input->candidates_total);
  if (edit_distance <= filter_input->max_error) {
    ++(filter_input->candidates_tp);
  } else {
    ++(filter_input->candidates_tn);
  }
  // Free
  edit_bpm_pattern_free(&bpm_pattern,filter_input->mm_allocator);
}


