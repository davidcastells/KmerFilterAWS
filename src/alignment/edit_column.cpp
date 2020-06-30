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
 * DESCRIPTION: Dynamic-programming alignment table
 */

#include "edit_column.h"

/*
 * DP-Table Setup
 */
void edit_column_allocate(
    edit_column_t* const edit_column,
    const int pattern_length,
    const int text_length,
    mm_allocator_t* const mm_allocator) {
  // Allocate DP table
  edit_column->num_rows = pattern_length + 1;
  const int num_columns = text_length + 1;
  edit_column->num_columns = num_columns;
  edit_column->column_curr = mm_allocator_calloc(mm_allocator,pattern_length+1,int,false);
  edit_column->column_prev = mm_allocator_calloc(mm_allocator,pattern_length+1,int,false);
}
void edit_column_free(
    edit_column_t* const edit_column,
    mm_allocator_t* const mm_allocator) {
  mm_allocator_free(mm_allocator,edit_column->column_curr);
  mm_allocator_free(mm_allocator,edit_column->column_prev);
}

