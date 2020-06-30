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

#include "../filter/pattern_tiling.h"

/*
 * Debug
 */
#define DEBUG_PATTERN_TILE_POSITION

/*
 * Pattern tiling
 */
void pattern_tiling_init(
    pattern_tiling_t* const pattern_tiling,
    const uint64_t pattern_length,
    const uint64_t pattern_tile_length,
    const uint64_t sequence_length,
    const uint64_t max_error) {
  // Calculate default tile dimensions & position
  pattern_tiling->pattern_max_error = max_error;
  pattern_tiling->pattern_tile_tall = pattern_tile_length;
  pattern_tiling->pattern_remaining_length = pattern_length;
  const int64_t pattern_band_width = sequence_length - pattern_length + 2*max_error;
  if (pattern_band_width < 0) {
    // It's supposed to be an invalid case because the filtering-process should have
    // detected (sequence_length < pattern_length) and trim the pattern accordingly
    fprintf(stderr,"Pattern tiling went wrong: invalid dimensions");
    exit(1);
  }
  pattern_tiling->pattern_band_width = pattern_band_width;
  pattern_tiling->sequence_length = sequence_length;
  pattern_tiling->tile_offset = 0;
  // Calculate current tile dimensions (adjusted to initial conditions)
  pattern_tiling->tile_next_offset_inc = BOUNDED_SUBTRACTION(pattern_tile_length,max_error,0);
  pattern_tiling->tile_tall = MIN(pattern_tile_length,pattern_length);
  pattern_tiling->tile_wide = pattern_tiling->pattern_band_width + pattern_tiling->tile_next_offset_inc;
  if (pattern_tiling->tile_offset+pattern_tiling->tile_wide > pattern_tiling->sequence_length) {
    pattern_tiling->tile_wide = pattern_tiling->sequence_length - pattern_tiling->tile_offset;
  }
}
void pattern_tiling_next(
    pattern_tiling_t* const pattern_tiling) {
  // DEBUG
#ifdef DEBUG_PATTERN_TILE_POSITION
  fprintf(stderr,">> Tile (pos=%lu,len=%lu,tall=%lu)\n",
      pattern_tiling->tile_offset,pattern_tiling->tile_wide,
      pattern_tiling->tile_tall);
#endif
  // Update tile dimensions
  pattern_tiling->tile_offset += pattern_tiling->tile_next_offset_inc;
  pattern_tiling->pattern_remaining_length -= pattern_tiling->tile_tall;
  // Calculate current tile dimensions
  pattern_tiling->tile_next_offset_inc = pattern_tiling->tile_tall;
  pattern_tiling->tile_tall = MIN(pattern_tiling->tile_tall,pattern_tiling->pattern_remaining_length);
  pattern_tiling->tile_wide = pattern_tiling->pattern_band_width + pattern_tiling->tile_next_offset_inc;
  if (pattern_tiling->tile_offset+pattern_tiling->tile_wide > pattern_tiling->sequence_length) {
    pattern_tiling->tile_wide = pattern_tiling->sequence_length - pattern_tiling->tile_offset;
  }
}
//uint64_t pattern_tiling_bound_matching_path(pattern_tiling_t* const pattern_tiling) {
//  if (pattern_tiling->prev_tile_match_position!=UINT64_MAX) {
//    const int64_t prev_tile_match_position = pattern_tiling->prev_tile_match_position;
//    const int64_t tile_match_position = pattern_tiling->tile_match_column + pattern_tiling->tile_offset;
//    const int64_t tile_match_position_proyection = tile_match_position - pattern_tiling->tile_tall;
//    // Calculate plausible match-band limits
//    const int64_t tile_match_band_begin =
//        BOUNDED_SUBTRACTION(tile_match_position_proyection,pattern_tiling->tile_distance,0);
//    const int64_t tile_match_band_end =
//        BOUNDED_ADDITION(tile_match_position_proyection,pattern_tiling->tile_distance,pattern_tiling->sequence_length-1);
//    // Calculate differences
//    const int64_t band_begin_difference = ABS(prev_tile_match_position - tile_match_band_begin);
//    const int64_t band_end_difference = ABS(prev_tile_match_position - tile_match_band_end);
//    // Keep matching column
//    pattern_tiling->prev_tile_match_position = tile_match_position;
//    // Return bound
//    return MAX(band_begin_difference,band_end_difference);
//  } else {
//    // Keep matching column
//    pattern_tiling->prev_tile_match_position = pattern_tiling->tile_match_column;
//    return 0;
//  }
//}
