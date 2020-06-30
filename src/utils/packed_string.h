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
 * DESCRIPTION: Encoded string in 3-bits bitmap (packed)
 */

#ifndef PACKED_STRING_H_
#define PACKED_STRING_H_

/*
 * Includes
 */
#include "../utils/commons.h"
#include "../system/mm_allocator.h"

/*
 * Masks
 */
// 1101 1011 0110 1101 1011 0110 1101 1011 0110 1101 1011 0110 1101 1011 0110 1101
//    D    B    6   D     B    6    D    B    6    D    B    6    D    B    6    D
#define PACKED_STRING_PADDING_BLOCK_X  0xDB6DB6DB6DB6DB6Dull
// 1110 1101 1011 0110 1101 1011 0110 1101 1011 0110 1101 1011 0110 1101 1011 0110
//    E    D    B    6    D    B    6    D    B    6    D    B    6    D    B    6
#define PACKED_STRING_PADDING_BLOCK_Y  0xEDB6DB6DB6DB6DB6ull
#define PACKED_STRING_LEADING_ONE      0x8000000000000000ull

/*
 * Dimensions
 */
#define NUM_BITS_PER_LETTER            3
#define NUM_LETTERS_PER_64WORD         (UINT64_LENGTH/NUM_BITS_PER_LETTER)
#define BITMAT_NUM_BLOCKS              NUM_LETTERS_PER_64WORD

/*
 * Packed String
 */
typedef struct {
  void* memory[BITMAT_NUM_BLOCKS];
  uint64_t* blocks_offseted[BITMAT_NUM_BLOCKS];
  // MM
  mm_allocator_t* mm_allocator;
} packed_string_t;

/*
 * Setup
 */
packed_string_t* packed_string_new(
    const char* const buffer,
    const int buffer_length,
    const int begin_padding_length,
    const int end_padding_length,
    const uint64_t padding_block,
    mm_allocator_t* const mm_allocator);
void packed_string_delete(
    packed_string_t* const packed_string);


#endif /* PACKED_STRING_H_ */
