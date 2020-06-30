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
 * DESCRIPTION: Simple Hash Implementation (String Key, Generic Value)
 */

#ifndef SHASH_H_
#define SHASH_H_

#include "../utils/commons.h"
#include "../system/mm_allocator.h"
#include "../utils/uthash.h"

/*
 * String Key Hash
 */
typedef struct {
  char* key;
  void* element;
  UT_hash_handle hh;
} shash_element_t;
typedef struct {
  shash_element_t* head;
  mm_allocator_t* mm_allocator;
} shash_t;
typedef struct {
  shash_t* shash;
  shash_element_t* current;
  shash_element_t* next;
} shash_iterator_t;

/*
 * Constructor
 */
shash_t* shash_new(mm_allocator_t* const mm_allocator);
void shash_clear(shash_t* const shash);
void shash_delete(shash_t* const shash);

/*
 * Handlers (Type-unsafe Accessors)
 */
shash_element_t* shash_handler_get_element(
    shash_t* const shash,
    char* const key);
void* shash_handler_get(
    shash_t* const shash,
    char* const key);
void shash_handler_insert(
    shash_t* const shash,
    char* const key,
    const int key_length,
    void* const element);
void shash_handler_remove(
    shash_t* const shash,
    char* const key);

/*
 * Accessors
 *   The @element is inserted by reference (not copied, nor freed at removal)
 *   In case of duplicates, the last element added remains
 */
#define shash_get(shash,key,type) ((type*)shash_handler_get(shash,key))
#define shash_insert(shash,key,key_length,element) shash_handler_insert(shash,key,key_length,(void*)element)
#define shash_remove(shash,key) shash_handler_remove(shash,key)

bool shash_is_contained(
    shash_t* const shash,
    char* const key);
int shash_get_num_elements(shash_t* const shash);

/*
 * Iterator
 */
#define SHASH_BEGIN_ITERATE(shash,it_skey,it_element,type) { \
  shash_element *shash_##sh_element, *shash_##tmp; \
  HASH_ITER(hh,shash,shash_##sh_element,shash_##tmp) { \
    type* const it_element = (type*)(shash_##sh_element->element); \
    char* const it_skey = shash_##sh_element->key;

#define SHASH_BEGIN_ELEMENT_ITERATE(shash,it_element,type) { \
  shash_element *shash_##sh_element, *shash_##tmp; \
  HASH_ITER(hh,shash,shash_##sh_element,shash_##tmp) { \
    type* const it_element = (type*)(shash_##sh_element->element);

#define SHASH_BEGIN_KEY_ITERATE(shash,it_skey) { \
  shash_element *shash_##sh_element, *shash_##tmp; \
  HASH_ITER(hh,shash,shash_##sh_element,shash_##tmp) { \
    char* const it_skey = shash_##sh_element->key;

#define SHASH_END_ITERATE }}

shash_iterator_t* shash_iterator_new(shash_t* const shash);
void shash_iterator_delete(shash_iterator_t* const iterator);

bool shash_iterator_eoi(shash_iterator_t* const iterator);
bool shash_iterator_next(shash_iterator_t* const iterator);
char* shash_iterator_get_key(shash_iterator_t* const iterator);
void* shash_iterator_get_element(shash_iterator_t* const iterator);

#endif /* SHASH_H_ */
