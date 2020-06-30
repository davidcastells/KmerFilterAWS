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
 * PROJECT: Wavefront Alignments Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Dynamic-programming alignment algorithm for Levenshtein distance (edit)
 */

#include "edit_dp.h"

/*
 * Edit distance computation [Ends-free] [DP-Columns]
 */
int edit_dp_distance(
    edit_column_t* const edit_column,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length) {
  // Parameters
  int* column_curr = edit_column->column_curr;
  int* column_prev = edit_column->column_prev;
  int h, v, min_distance = INT_MAX;
  // Initialize
  column_curr[0] = 0; // [Ends-free]
  for (v=0;v<=pattern_length;++v) column_prev[v] = v;
  // Compute DP
  for (h=1;h<=text_length;++h) {
    for (v=1;v<=pattern_length;++v) {
      int min = column_prev[v-1] + (text[h-1]!=pattern[v-1]);
      min = MIN(min,column_prev[v]+1); // Ins
      min = MIN(min,column_curr[v-1]+1); // Del
      column_curr[v] = min;
    }
    // Check min distance
    const int new_min_distance = column_curr[pattern_length];
    if (new_min_distance < min_distance) min_distance = new_min_distance;
    // Swap
    SWAP(column_curr,column_prev);
  }
  // Return distance
  return min_distance;
}
/*
 * Edit distance computation using raw DP-Table (banded)
 */
int edit_dp_distance_banded(
    edit_column_t* const edit_column,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int bandwidth) {
  // Parameters
  int* column_curr = edit_column->column_curr;
  int* column_prev = edit_column->column_prev;
  int h, v, min_distance = INT_MAX;
  // Compute band limits (account for mandatory indels)
  int text_band_init = bandwidth + ((text_length>pattern_length) ? (text_length-pattern_length) : 0);
  if (text_band_init > text_length) text_band_init = text_length;
  int pattern_band_init = bandwidth + ((pattern_length>text_length) ? (pattern_length-text_length) : 0);
  if (pattern_band_init > pattern_length) pattern_band_init = pattern_length;
  // Initialize [Text Ends-free]
  column_curr[0] = 0;
  for (v=0;v<=pattern_band_init;++v) column_prev[v] = v;
  // Compute DP
  int lo_band = 1;
  int hi_band = pattern_band_init;
  for (h=1;h<=text_length;++h) {
    // Compute band limits
    if (h > text_band_init) {
      column_curr[lo_band] = INT16_MAX; // Init corner case
      ++lo_band;
    }
    if (hi_band < pattern_length) {
      ++hi_band;
      column_prev[hi_band] = INT16_MAX; // Init corner case
    }
    // Compute column
    for (v=lo_band;v<=hi_band;++v) {
      int min = column_prev[v-1] + (text[h-1]!=pattern[v-1]);
      min = MIN(min,column_prev[v]+1); // Ins
      min = MIN(min,column_curr[v-1]+1); // Del
      column_curr[v] = min;
    }
    // Check min distance
    if (hi_band == pattern_length) {
      const int new_min_distance = column_curr[pattern_length];
      if (new_min_distance < min_distance) min_distance = new_min_distance;
    }
    // Swap
    SWAP(column_curr,column_prev);
  }
  // Return distance
  return min_distance;
}
