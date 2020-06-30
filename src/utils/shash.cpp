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

#include "../utils/shash.h"

/*
 * Allocators
 */
#undef  uthash_malloc
#define uthash_malloc(sz) (shash->mm_allocator!=NULL ? mm_allocator_malloc(shash->mm_allocator,sz) : malloc(sz))
#undef  uthash_free
#define uthash_free(ptr,sz) (shash->mm_allocator!=NULL) ? mm_allocator_free(shash->mm_allocator,ptr) : free(ptr);

/*
 * Constructor
 */
shash_t* shash_new(mm_allocator_t* const mm_allocator) {
  shash_t* const shash = (mm_allocator==NULL) ?
      (shash_t*) malloc(sizeof(shash_t)) :
      mm_allocator_alloc(mm_allocator,shash_t);
  shash->head = NULL; // uthash initializer
  shash->mm_allocator = mm_allocator;
  return shash;
}
void shash_clear(shash_t* shash) {
  shash_element_t *shash_element, *tmp;
  HASH_ITER(hh,shash->head,shash_element,tmp) {
    HASH_DEL(shash->head,shash_element);
    uthash_free(shash_element->key,0);
    uthash_free(shash_element,0);
  }
  shash->head = NULL; // uthash initializer
}
void shash_delete(shash_t* shash) {
  shash_clear(shash);
  uthash_free(shash,0);
}
/*
 * Handlers (Type-unsafe Accessors)
 */
shash_element_t* shash_handler_get_element(
    shash_t* const shash,
    char* const key) {
  shash_element_t* shash_element;
  HASH_FIND_STR(shash->head,key,shash_element);
  return shash_element;
}
void* shash_handler_get(
    shash_t* const shash,
    char* const key) {
  shash_element_t* shash_element = shash_handler_get_element(shash,key);
  return (shash_element!=NULL) ? shash_element->element : NULL;
}
void shash_handler_insert(
    shash_t* const shash,
    char* const key,
    const int key_length,
    void* const element) {
  shash_element_t* shash_element = shash_handler_get_element(shash,key);
  if (shash_element==NULL) {
    // Allocate
    shash_element = (shash_element_t*) uthash_malloc(sizeof(shash_element_t));
    // Set key/element
    shash_element->key = (char*) uthash_malloc(key_length+1);
    memcpy(shash_element->key,key,key_length);
    shash_element->key[key_length+1] = '\0';
    shash_element->element = element;
    // Add to hash
    HASH_ADD_KEYPTR(hh,shash->head,shash_element->key,key_length,shash_element);
  } else {
    // Replace element (no free whatsoever)
    shash_element->element = element;
  }
}
void shash_handler_remove(
    shash_t* shash,
    char* const key) {
  shash_element_t* shash_element = shash_handler_get_element(shash,key);
  if (shash_element) {
    HASH_DEL(shash->head,shash_element);
    uthash_free(shash_element->key,0);
    uthash_free(shash_element,0);
  }
}
/*
 * Accessors
 */
bool shash_is_contained(
    shash_t* const shash,
    char* const key) {
  return (shash_handler_get_element(shash,key)!=NULL);
}
int shash_get_num_elements(shash_t* const shash) {
  return (int)HASH_COUNT(shash->head);
}
/*
 * Iterator
 */
shash_iterator_t* shash_iterator_new(shash_t* const shash) {
  // Allocate
  shash_iterator_t* const iterator = (shash_iterator_t*) uthash_malloc(sizeof(shash_iterator_t));
  // Init
  iterator->shash = shash;
  iterator->current = NULL;
  iterator->next = shash->head;
  return iterator;
}
void shash_iterator_delete(shash_iterator_t* const iterator) {
  shash_t* const shash = iterator->shash;
  uthash_free(iterator,0);
}
bool shash_iterator_eoi(shash_iterator_t* const iterator) {
  return (iterator->next != NULL);
}
bool shash_iterator_next(shash_iterator_t* const iterator) {
  if (iterator->next != NULL) {
    iterator->current = iterator->next;
    iterator->next = (shash_element_t*) iterator->next->hh.next;
    return true;
  } else {
    iterator->current = NULL;
    return false;
  }
}
char* shash_iterator_get_key(shash_iterator_t* const iterator) {
  return iterator->current->key;
}
void* shash_iterator_get_element(shash_iterator_t* const iterator) {
  return iterator->current->element;
}

