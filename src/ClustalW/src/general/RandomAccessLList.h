/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * Author: Mark Larkin UCD Conway Institute
 *
 * This class is to be used as a linked list type container. 
 * The features are as follows:
 * 1) It allows efficient deletions anywhere.
 * 2) It also allows efficient random access.
 * 3) Its maximum size is fixed at creation.
 * 4) Items can only be added at the end.
 * 5) when it is destroyed it will delete the structure, but not the elements.
 *
 * Limitations:
 * It returns the ListElements directly, this is not the best way. A better way would be
 * to return an iterator of some kind, but this would take alot more work. At the moment
 * it is possible to have 2 LList objects and get an element from one list, and delete
 * it from another list. But this is ok for the time being. I will fix this later.
 */
#ifndef RANDOMACCESSLLIST_H
#define RANDOMACCESSLLIST_H

namespace clustalw
{

template <class T>
class ListElement
{
    public:
        ListElement()
        {
            // Do nothing
        }
        
        ListElement(T* e, ListElement* n, ListElement* p, unsigned int _index)
        {
            element = e;
            next = n;
            prev = p;
            index = _index;
        }
                
        ~ListElement()
        {
            element = 0;
            next = 0;
        }
        
        unsigned int getIndex(){return index;}
                
        T* element;
        ListElement* next;
        ListElement* prev;
    private:    
        unsigned int index;
};

template <class T>
class RandomAccessLList
{
    public:
        RandomAccessLList(int size)
         : maxSize(size),
           numElements(0)
        {
            elements = new ListElement<T>*[size];
            for(int i = 0; i < size; i++)
            {
                elements[i] = 0;
            }
        }
        
        ~RandomAccessLList()
        {
            for(int i = 0; i < maxSize; i++)
            {
                if(elements[i] != 0)
                {
                    delete elements[i];
                }
            }
            delete [] elements;        
        }
        
        unsigned int getNumElements(){return numElements;}
        
        void addElementToBack(T* element)
        {
            if(numElements < maxSize)
            {
                // Add the element
                elements[numElements] = new ListElement<T>(element, 0, 0, numElements);
                
                if(numElements == 0)
                {
                    // First element
                    elements[numElements]->next = 0;
                    elements[numElements]->prev = 0;
                    firstElementInList = elements[numElements];
                    numElements++;             
                }
                else if(numElements > 0) // If it is not the first element
                {
                    elements[numElements - 1]->next = elements[numElements];
                    elements[numElements]->prev = elements[numElements - 1];
                    elements[numElements]->next = 0;
                    numElements++;
                }
            }
        }
        
        ListElement<T>* getAt(int index)
        {
            if(index >= 0 && index < maxSize)
            {
                return elements[index];
            }
            else
            {
                cout << "\nCannot get item : " << index << "\n";
                return 0;
            }
        }
        
        ListElement<T>* getFirst()
        {
            return firstElementInList;
        }        
        
        void removeItem(ListElement<T>* item)
        {
            if(item != 0 )
            {
                int index = item->getIndex();
                if(elements[index] != 0)
                { 
                    if(item->prev == 0 && item->next == 0) // Last item in the list
                    {
                        firstElementInList = 0;
                        delete item; // Not sure if this is good or not.
                        item = 0;
                        numElements--;
                    }
                    else if(item->prev == 0) // remove First in the list
                    {
                        firstElementInList = item->next;
                        firstElementInList->prev = 0;
                        item->next = 0;
                        delete item; // Not sure if this is good or not.
                        item = 0;
                        numElements--;
                    }
                    else if(item->next == 0) // remove last item
                    {
                        ListElement<T>* prevElem = item->prev;
                        prevElem->next = 0;
                        item->prev = 0;
                        delete item; // Not sure if this is good or not.
                        item = 0;
                        numElements--;                
                    }
                    else // remove middle item
                    {
                        ListElement<T>* prevElem = item->prev;
                        ListElement<T>* nextElem = item->next;
                        prevElem->next = nextElem;
                        nextElem->prev = prevElem;
                        item->prev = 0;
                        item->next = 0;
                        delete item; // Not sure if this is good or not.
                        item = 0;
                        numElements--;                 
                    }
                    elements[index] = 0;
                }
            }
            else
            {
                cout << "ERROR: trying to remove item outside the bounds\n";
            }            
        }
        
    private:    
        ListElement<T>** elements;
        ListElement<T>* firstElementInList;
        int maxSize;
        unsigned int numElements;
        
};

}

#endif
