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

#include "../utils/packed_string.h"
#include "../system/mm_allocator.h"

/*
 * Setup
 */
packed_string_t* packed_string_new(
    const char* const buffer,
    const int buffer_length,
    const int begin_padding_length,
    const int end_padding_length,
    const uint64_t padding_block,
    mm_allocator_t* const mm_allocator) {
  // Allocate
  packed_string_t* const packed_string = mm_allocator_alloc(mm_allocator,packed_string_t);
  packed_string->mm_allocator = mm_allocator;
  // Compute sizes
  const int begin_padding_64words = DIV_CEIL(begin_padding_length*NUM_BITS_PER_LETTER,NUM_LETTERS_PER_64WORD);
  const int buffer_64words = DIV_CEIL(buffer_length*NUM_BITS_PER_LETTER,NUM_LETTERS_PER_64WORD);
  const int end_padding_64words = DIV_CEIL(end_padding_length*NUM_BITS_PER_LETTER,NUM_LETTERS_PER_64WORD);
  // Allocate buffers
  const int string_enc_size = (begin_padding_64words + buffer_64words + end_padding_64words)*UINT64_SIZE;
  int offset;
  for (offset=0;offset<BITMAT_NUM_BLOCKS;++offset) {
    packed_string->memory[offset] = mm_allocator_calloc(mm_allocator,string_enc_size,char,true);
    packed_string->blocks_offseted[offset] = ((uint64_t*)packed_string->memory[offset]) + begin_padding_64words;
  }
  // Add begin padding
  int i;
  for (offset=0;offset<BITMAT_NUM_BLOCKS;++offset) {
    uint64_t* const string_mem = (uint64_t*) packed_string->memory[offset];
    for (i=0;i<begin_padding_64words;++i) {
      string_mem[i] = padding_block;
    }
  }
  // Encode buffer (offset 0)
  uint64_t* blocks_offseted_0 = packed_string->blocks_offseted[0];
  uint64_t block = 0;
  for (i=0,offset=0;i<buffer_length;++i) {
    // Encode letter (from right to left)
    block = block | (((uint64_t)buffer[i]) << offset);
    // Increment offset
    offset += NUM_BITS_PER_LETTER;
    if (offset == 63) {
      // Write & Next
      *blocks_offseted_0 = PACKED_STRING_LEADING_ONE | block;
      ++blocks_offseted_0;
      // Reset
      block = 0;
      offset = 0;
    }
  }
  packed_string->blocks_offseted[0][buffer_64words] = padding_block; // Pad last
  // Shift buffer (offset 1-20)
  int offset_block, offset_next_block;
  uint64_t next_block = 0;
  for (i=0;i<buffer_64words;++i) {
    // Fetch block and next block
    block = packed_string->blocks_offseted[0][i];
    next_block = packed_string->blocks_offseted[0][i+1];
    offset_block = 0;
    offset_next_block = 63;
    // Shift
    for (offset=0;offset<BITMAT_NUM_BLOCKS;++offset) {
      packed_string->blocks_offseted[offset][i] =
          PACKED_STRING_LEADING_ONE | ((block >> offset_block) | (next_block << offset_next_block));
      offset_block += 3;
      offset_next_block -= 3;
    }
  }
  // Add end padding
  for (offset=0;offset<8;++offset) {
    uint64_t* const string_mem = (uint64_t*) (packed_string->memory[offset] + begin_padding_64words + buffer_64words);
    for (i=0;i<end_padding_64words;++i) {
      string_mem[i] = padding_block;
    }
  }
  // Return
  return packed_string;
}
void packed_string_delete(
    packed_string_t* const packed_string) {
  int offset;
  for (offset=0;offset<8;++offset) {
    mm_allocator_free(packed_string->mm_allocator,packed_string->memory[offset]);
  }
  mm_allocator_free(packed_string->mm_allocator,packed_string);
}


