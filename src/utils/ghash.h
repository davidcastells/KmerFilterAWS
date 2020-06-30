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
 * DESCRIPTION: Simple Hash Implementation (Generic Key, Generic Value)
 */

#ifndef GHASH_H_
#define GHASH_H_

#include "../utils/commons.h"
#include "../system/mm_allocator.h"
#include "../utils/uthash.h"

/*
 * General Key Hash
 */
typedef struct {
  void* key;
  void* element;
  UT_hash_handle hh;
} ghash_element_t;
typedef struct {
  ghash_element_t* head;
  int key_size;
  mm_allocator_t* mm_allocator;
} ghash_t;
typedef struct {
  ghash_t* ghash;
  ghash_element_t* current;
  ghash_element_t* next;
} ghash_iterator_t;

/*
 * Constructor
 */
ghash_t* ghash_new(
    const int key_size,
    mm_allocator_t* const mm_allocator);
void ghash_clear(ghash_t* const ghash);
void ghash_delete(ghash_t* const ghash);

/*
 * Handlers (Type-unsafe Accessors)
 */
ghash_element_t* ghash_handler_get_element(
    ghash_t* const ghash,
    const void* const key);
void* ghash_handler_get(
    ghash_t* const ghash,
    const void* const key);
void ghash_handler_insert(
    ghash_t* const ghash,
    const void* const key,
    void* const element);
void ghash_handler_remove(
    ghash_t* const ghash,
    const void* const key);

/*
 * Accessors
 *   The @element is inserted by reference (not copied, nor freed at removal)
 *   In case of duplicates, the last element added remains
 */
#define ghash_get(ghash,key,type) ((type*)ghash_handler_get(ghash,(void*)(key)))
#define ghash_insert(ghash,key,element) ghash_handler_insert(ghash,(void*)(key),(void*)(element))
#define ghash_remove(ghash,key) ghash_handler_remove(ghash,(void*)(key))

bool ghash_is_contained(
    ghash_t* const ghash,
    const void* const key);
int ghash_get_num_elements(ghash_t* const ghash);
int ghash_get_size(ghash_t* const ghash);

/*
 * Iterator
 */
#define GHASH_BEGIN_ITERATE(ghash,it_gkey,it_element,type) { \
  ghash_element_t *ghash_##ih_element, *ghash_##tmp; \
  HASH_ITER(hh,ghash->head,ghash_##ih_element,ghash_##tmp) { \
    type* const it_element = (type*)(ghash_##ih_element->element);
#define GHASH_END_ITERATE }}

#define GHASH_BEGIN_ITERATE__KEY(ghash,it_gkey,it_element,type) { \
  ghash_element_t *ghash_##ih_element, *ghash_##tmp; \
  HASH_ITER(hh,ghash->head,ghash_##ih_element,ghash_##tmp) { \
    type* const it_element = (type*)(ghash_##ih_element->element); \
    void* const it_gkey = ghash_##ih_element->key;
#define GHASH_END_ITERATE__KEY }}

ghash_iterator_t* ghash_iterator_new(ghash_t* const ghash);
void ghash_iterator_delete(ghash_iterator_t* const ghash_iterator);

bool ghash_iterator_eoi(ghash_iterator_t* const iterator);
bool ghash_iterator_next(ghash_iterator_t* const ghash_iterator);
void* ghash_iterator_get_key(ghash_iterator_t* const ghash_iterator);
void* ghash_iterator_get_element(ghash_iterator_t* const ghash_iterator);

#endif /* GHASH_H_ */
