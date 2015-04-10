/* -*- mode: c; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */

/*********************************************************************
 * Clustal Omega - Multiple sequence alignment
 *
 * Copyright (C) 2010 University College Dublin
 *
 * Clustal-Omega is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This file is part of Clustal-Omega.
 *
 ********************************************************************/

/*
 * RCS $Id: util-C.h 155 2010-11-17 12:18:47Z fabian $
 */

// Utility subroutines


#ifndef MAIN
#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream 
#include <cstdio>     // printf
#include <stdlib.h>   // exit
#include <time.h>     // clock
#endif
#include <sys/time.h>

/////////////////////////////////////////////////////////////////////////////////////
// Arithmetics
/////////////////////////////////////////////////////////////////////////////////////

//// max and min
inline double dmax(double x, double y) { return (x>y? x : y);}
inline double dmin(double x, double y) { return (x<y? x : y);}
inline int imax(int x, int y) { return (x>y? x : y);}
inline int imin(int x, int y) { return (x<y? x : y);}
inline int iabs(int x) { return (x>=0? x : -x);}

// Rounding up, rounding down and rounding to nearest integer
inline int iceil(double x)  {return int(ceil(x));}
inline int ifloor(double x) {return int(floor(x));}
inline int iround(double x) {return int(floor(x+0.5));}

//// Generalized mean: d=0: sqrt(x*y)  d=1: (x+y)/2  d->-inf: min(x,y)  d->+inf: max(x,y)  
inline double fmean(double x, double y, double d) { return pow( (pow(x,d)+pow(y,d))/2 ,1./d);}

// log base 2
inline float log2(float x)  {return (x<=0? (float)(-100000):1.442695041*log(x));}
inline float log10(float x) {return (x<=0? (float)(-100000):0.434294481*log(x));}


/////////////////////////////////////////////////////////////////////////////////////
// fast log base 2
/////////////////////////////////////////////////////////////////////////////////////

// This function returns log2 with a max abolute deviation of +/- 1.5E-5 (typically 0.8E-5). 
// It takes 1.42E-8 s  whereas log2(x) takes 9.5E-7 s. It is hence 9.4 times faster. 
// It makes use of the representation of 4-byte floating point numbers:
// seee eeee emmm mmmm mmmm mmmm mmmm mmmm
// s is the sign, 
// the following 8 bits, eee eee e, give the exponent + 127 (in hex: 0x7f).
// The following 23 bits, m, give the mantisse, the binary digits behind the decimal point.
// In summary: x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeee-127) 
// The expression (((*(int *)&x) & 0x7f800000 ) >>23 )-0x7f is the exponent eeeeeeee, i.e. 
// the largest integer that is smaller than log2(x) (e.g. -1 for 0.9). *(int *)&x is an integer which 
// contains the bytes as the floating point variable x is represented in memory.
// Check:  assert( sizeof(f) == sizeof(int) );
// Check:  assert( sizeof(f) == 4 );
inline float fast_log2(float x) 
{
  static float lg2[1025];         // lg2[i] = log2[1+x/1024]
  static float diff[1025];        // diff[i]= (lg2[i+1]-lg2[i])/8096 (for interpolation)
  static char initialized;
  if (x<=0) return -100000;
  if (!initialized)   //First fill in the arrays lg2[i] and diff[i]
    {
      float prev = 0.0f;
      lg2[0] = 0.0f;
      for (int i=1; i<=1024; ++i) 
	{
	  lg2[i] = log(float(1024+i))*1.442695041-10.0f;
	  diff[i-1] = (lg2[i]-prev)*1.2352E-4;
	  prev = lg2[i];
	}
      initialized=1;
    }  
  int a = (((*((int *)&x)) & 0x7F800000) >>23 )-0x7f;
  int b =  ((*((int *)&x)) & 0x007FE000) >>13;
  int c =  ((*((int *)&x)) & 0x00001FFF);
  return a + lg2[b] + diff[b]*(float)(c);
}

/////////////////////////////////////////////////////////////////////////////////////
// fast 2^x
// ATTENTION: need to compile with g++ -fno-strict-aliasing when using -O2 or -O3!!!
// Relative deviation < 1.5E-4
/////////////////////////////////////////////////////////////////////////////////////
inline float fpow2(float x)
{
  if (x>=128) return FLT_MAX;
  if (x<=-128) return FLT_MIN;
  int *px = (int*)(&x);                 // store address of float as pointer to long
  float tx = (x-0.5f) + (3<<22);        // temporary value for truncation: x-0.5 is added to a large integer (3<<22)
  int lx = *((int*)&tx) - 0x4b400000;   // integer value of x
  float dx = x-(float)(lx);             // float remainder of x
  x = 1.0f + dx*(0.6960656421638072f          // cubic apporoximation of 2^x
	   + dx*(0.224494337302845f           // for x in the range [0, 1]
           + dx*(0.07944023841053369f)));
  *px += (lx<<23);                            // add integer power of 2 to exponent
  return x;
}

/////////////////////////////////////////////////////////////////////////////////////
// ATTENTION: 
// Can't be used with -O2/-O3 optimization on some compilers !
// Works with g++ version 4.1, but not with 3.4, in which case it returns values 
// that are a factor 1.002179942 too low 
//
// Fast pow2 routine (Johannes Soeding)  
// Same speed as fpow2(), but *relative* deviation < 1.2E-7
// Makes use of the binary representation of floats in memory:
//   x = (-1)^s * 1.mmmmmmmmmmmmmmmmmmmmmm * 2^(eeeeeee-127) 
// is represented as
// 31        23                   7654 3210
//  seee eeee emmm mmmm mmmm mmmm mmmm mmmm
// s is the sign, the 8 bits eee eee e are the exponent + 127 (in hex: 0x7f), 
// and the following 23 bits m give the mantisse.
// We decompose the argument x = a + b, with integer a and 0 <= b < 1
// Therefore 2^x = 2^a * 2^b  where a is the binary exponent of 2^x 
// and  1 <= 2^b < 2,  i.e. 2^b determines the mantisse uniquely. 
// To calculate 2^b, we split b into the first 10 bits and the last 13 bits,
// b = b' + c, and then look up the mantisse of 2^b' in a precomputed table.
// We use the residual c to interpolate between the mantisse for 2^b' and 2(b'+1/1024) 
/////////////////////////////////////////////////////////////////////////////////////
inline float fast_pow2(float x) 
{
  if (x<=-127) return 5.9E-39;
  if (x>=128)  return 3.4E38;
  static char initialized=0;
  static unsigned int pow2[1025];
  static unsigned int diff[1025];
  static int y = 0;
  if (!initialized)   //First fill in the pow2-vector
    {
      float f;
      unsigned int prev = 0;
      pow2[0] = 0;
      for (int b=1; b<1024; b++) 
	{
	  f=pow(2.0,float(b)/1024.0);  
	  pow2[b]=(*((unsigned int *)(&f)) & 0x7FFFFF); // store the mantisse of 2^(1+b/1024)
	  diff[b-1]=pow2[b]-prev;
	  prev=pow2[b];
	}
      pow2[1024]=0x7FFFFF;
      diff[1023]=pow2[1024]-prev;
      initialized=1;
    }  
      
  int *px = (int *)(&x);                              // store address of float as pointer to int
  int E = ((*px & 0x7F800000)>>23)-127;               // E is exponent of x and is <=6
  unsigned int M=(*px & 0x007FFFFF) | 0x00800000;     // M is the mantisse 1.mmm mmmm mmmm mmmm mmmm mmmm
  int a,b,c;
  if (x>=0) 
    {
      if (E>=0) {
	a = 0x3F800000 + ((M<<E) & 0x7F800000);       // a is exponent of 2^x, beginning at bit 23 
	b = ((M<<E) & 0x007FE000)>>13;                
	c = ((M<<E) & 0x00001FFF);
      } else {
	a =  0x3F800000;                              // a = exponent of 2^x = 0
	b = ((M>>(-E)) & 0x007FE000)>>13;          
	c = ((M>>(-E)) & 0x00001FFF);
      }
    } 
  else 
    {
      if (E>=0) {
	a = 0x3F000000 - ((M<<E) & 0x7F800000);       // a is exponent of 2^x
	b = (0x00800000-(int)((M<<E) & 0x007FFFFF)) >>13;          
	c = (0x00800000-(int)((M<<E) & 0x007FFFFF)) & 0x00001FFF;          
      } else {
	a = 0x3F000000;                               // a = exponent of 2^x = -1
	b = (0x00800000-(int)((M>>(-E)) & 0x007FFFFF)) >>13;          
	c = (0x00800000-(int)((M>>(-E)) & 0x007FFFFF)) & 0x00001FFF;          
      }
    }
/*   printf("x=%0X\n",*px); */
/*   printf("E=%0X\n",E); */
/*   printf("M=%0X\n",M); */
/*   printf("a=%0X\n",a); */
/*   printf("b=%0X\n",b); */
  y = a | (pow2[b] + ((diff[b]*c)>>13) );
  /*   printf("2^x=%0X\n",*px); */
  return *((float*)&y);
}



// Normalize a float array such that it sums to one
// If it sums to 0 then assign def_array elements to array (optional)
inline float NormalizeTo1(float* array, int length, float* def_array=NULL) 
{
  float sum=0.0f;
  int k;
  for (k=0; k<length; k++) sum+=array[k];
  if (sum!=0.0f) 
    {
      float fac=1.0/sum;
      for (k=0; k<length; k++) array[k]*=fac;
    } 
  else if (def_array)
    for (k=0; k<length; k++) array[k]=def_array[k];
  return sum;
}

// Normalize a float array such that it sums to x
// If it sums to 0 then assign def_array elements to array (optional)
inline float NormalizeToX(float* array, int length, float x, float* def_array=NULL) 
{
  float sum=0.0;
  int k;
  for (k=0; k<length; k++) sum+=array[k];
  if (sum) 
    {
      float fac=x/sum;
      for (k=0; k<length; k++) array[k]*=fac;
    } 
  else if (def_array)
    for (k=0; k<length; k++) array[k]=def_array[k];
  return sum;
}

/////////////////////////////////////////////////////////////////////////////////////
// Similar to spintf("%*g,w,val), but displays maximum number of digits within width w
/////////////////////////////////////////////////////////////////////////////////////
inline char* sprintg(float val, int w)
{
  static char str[100];
  float log10val = log10(fabs(val));
  int neg = (val<0? 1: 0);
  if (log10val >= w-neg-1 || -log10val > 3) 
    {
      // positive exponential 1.234E+06
      // negative exponential 1.234E-06
      int d = w-6-neg;
      sprintf(str,"%*.*e",w,d<1?1:d,val);
    }
  else 
    {
      int d = log10val>0? w-2-neg-int(log10val): w-2-neg;
      sprintf(str,"%#*.*f",w,d,val);
    }
  return str;
}

/////////////////////////////////////////////////////////////////////////////////////
// String utilities
/////////////////////////////////////////////////////////////////////////////////////

//the integer. If no integer is found, returns INT_MIN and sets ptr to NULL      /* MR1 */
inline int strtoi(const char*& ptr)
{
    int i;
    const char* ptr0=ptr;
    if (!ptr) return INT_MIN;
    while (*ptr!='\0' && !(*ptr>='0' && *ptr<='9')) ptr++;
    if (*ptr=='\0') {
        ptr=0;
        return INT_MIN;
    }
    if (*(ptr-1)=='-' && ptr>ptr0) i=-atoi(ptr); else i=atoi(ptr);
    while (*ptr>='0' && *ptr<='9') ptr++;
    return i;
}


//Same as strint, but interpretes '*' as default /* MR1 */
inline int strtoi_(const char*& ptr, int deflt=INT_MAX)
{
    int i;
    if (!ptr) return INT_MIN;
    while (*ptr!='\0' && !(*ptr>='0' && *ptr<='9') && *ptr!='*') ptr++;
    if (*ptr=='\0') {
        ptr=0;
        return INT_MIN;
    }
    if (*ptr=='*') {
        ptr++;
        return deflt;
    }
    if (*(ptr-1)=='-') i=atoi(ptr-1);
    else i=atoi(ptr);
    while (*ptr>='0' && *ptr<='9') ptr++;
    return i;
}


// Returns leftmost integer in ptr and sets the pointer to first char after 
// the integer. If no integer is found, returns INT_MIN and sets pt to NULL
int strint(char*& ptr)
{
  int i;
  char* ptr0=ptr;
  if (!ptr) return INT_MIN;
  while (*ptr!='\0' && !(*ptr>='0' && *ptr<='9')) ptr++;
  if (*ptr=='\0') 
    {
      ptr=0;
      return INT_MIN;
    }
  if (*(ptr-1)=='-' && ptr>ptr0) i=-atoi(ptr); else i=atoi(ptr);
  while (*ptr>='0' && *ptr<='9') ptr++;
  return i;
}

// Same as strint, but interpretes '*' as default
int strinta(char*& ptr, int deflt=99999)
{
  int i;
  if (!ptr) return INT_MIN;
  while (*ptr!='\0' && !(*ptr>='0' && *ptr<='9') && *ptr!='*') ptr++;
  if (*ptr=='\0') 
    {
      ptr=0;
      return INT_MIN;
    }
  if (*ptr=='*') 
    {
      ptr++;
      return deflt;
    }
  if (*(ptr-1)=='-') i=atoi(ptr-1);
  else i=atoi(ptr);
  while (*ptr>='0' && *ptr<='9') ptr++;
  return i;
}

// Returns leftmost float in ptr and sets the pointer to first char after  
// the float. If no float is found, returns FLT_MIN and sets pt to NULL /* MR1 */
float strflt(char*& ptr)
{
    float i;
    char* ptr0=ptr;
    if (!ptr) return FLT_MIN;
    while (*ptr!='\0' && !(*ptr>='0' && *ptr<='9')) ptr++;
    if (*ptr=='\0')
        {
            ptr=0;
            return FLT_MIN;
        }
    if (ptr>ptr0 && *(ptr-1)=='-') i=-atof(ptr); else i=atof(ptr);
    while ((*ptr>='0' && *ptr<='9') || *ptr=='.') ptr++;
    return i;
}

// Same as strint, but interpretes '*' as default  /* MR1 */
float strflta(char*& ptr, float deflt=99999)
{
    float i;
    if (!ptr) return FLT_MIN;
    while (*ptr!='\0' && !(*ptr>='0' && *ptr<='9') && *ptr!='*') ptr++;
    if (*ptr=='\0')
        {
            ptr=0;
            return FLT_MIN;
        }
    if (*ptr=='*')
        {
            ptr++;
            return deflt;
        }
    if (*(ptr-1)=='-') i=-atof(ptr);
    else i=atof(ptr);
    while ((*ptr>='0' && *ptr<='9') || *ptr=='.') ptr++;
    return i;
}


// Removes the newline and other control characters at the end of a string (if present)
// and returns the new length of the string (-1 if str is NULL)
inline int chomp(char str[]) 
{
  if (!str) return -1;
  int l=0;
  for (l=strlen(str)-1; l>=0 && str[l]<32; l--);
  str[++l]='\0';
  return l;
}

// Emulates the ifstream::getline method; similar to fgets(str,maxlen,FILE*), 
// but removes the newline at the end and returns NULL if at end of file or read error
inline char* fgetline(char str[], const int maxlen, FILE* file) 
{
  if (!fgets(str,maxlen,file)) return NULL;
  if (chomp(str)+1>=maxlen)    // if line is cut after maxlen characters...
    while (fgetc(file)!='\n'); // ... read in rest of line  
  return(str);
}

// copies substring str[a,b] into substr and returns substr 
char *substr(char* substr, char* str, int a, int b)
{
  if (b<a) {int i=b; b=a; a=i;}
  if (b-a>1000) 
    {printf("Function substr: >1000 chars to copy. Exiting.\n"); throw 6;}
  char* dest=substr;
  char* source=str+a;
  char* send=str+b;
  while (*source!='\0' && source<=send) *(dest++) = *(source++); 
  *dest='\0';
  return substr;
}


// Returns pointer to first non-white-space character in str OR to NULL if none found
inline char* strscn(char* str)
{
  if (!str) return NULL;
  char* ptr=str;
  while (*ptr!='\0' && *ptr<=32) ptr++;
  return (*ptr=='\0')? NULL: ptr;
}

// Returns pointer to first white-space character in str OR to NULL if none found   /* MR1 */
inline char* strscn_ws(char* str)
{
    if (!str) return NULL;
    char* ptr=str;
    while (*ptr!='\0' && *ptr>32) ptr++;
    return (*ptr=='\0')? NULL: ptr;
}

//Returns pointer to first non-white-space character in str OR to NULL if none found  /* MR1 */
inline const char* strscn_c(const char* str)
{
    if (!str) return NULL;
    const char* ptr=str;
    while (*ptr!='\0' && isspace(*ptr)) ptr++;
    return (*ptr=='\0') ? NULL : ptr;
}

// Returns pointer to first  non-white-space character in str OR to end of string '\0' if none found
inline char* strscn_(char* str)
{
  if (!str) return NULL;
  char* ptr=str;
  while (*ptr!='\0' && *ptr<=32) ptr++;
  return ptr;
}

// Returns pointer to first non-c character in str OR to NULL if none found
inline char* strscn(char* str, const char c)
{
  if (!str) return NULL;
  char* ptr=str;
  while (*ptr!='\0' && *ptr==c) ptr++;
  return (*ptr=='\0')? NULL: ptr;
}

// Returns pointer to first  non-c character in str OR to end of string '\0' if none found
inline char* strscn_(char* str, const char c)
{
  if (!str) return NULL;
  char* ptr=str;
  while (*ptr!='\0' && *ptr==c) ptr++;
  return ptr;
}

// Cuts string at first white space character found by overwriting it with '\0'. 
// Returns pointer to next non-white-space char OR to NULL if no such char found 
inline char* strcut(char* str)
{
  if (!str) return NULL;
  char* ptr=str;
  while (*ptr!='\0' && *ptr>32) ptr++;
  if (*ptr=='\0') return NULL;
  *ptr='\0';
  ptr++;
  while (*ptr!='\0' && *ptr<=32) ptr++;
  return (*ptr=='\0')? NULL:ptr;
}

// Cuts string at first white space character found by overwriting it with '\0'. 
// Returns pointer to next non-white-space char OR to end of string '\0' if none found
inline char* strcut_(char* str)
{
  if (!str) return NULL;
  char* ptr=str;
  while (*ptr!='\0' && *ptr>32) ptr++;
  if (*ptr=='\0') return ptr;
  *ptr='\0';
  ptr++;
  while (*ptr!='\0' && *ptr<=32) ptr++;
  return ptr;
}

// Cuts string at first occurence of charcter c, by overwriting it with '\0'.
// Returns pointer to next char not equal c, OR to NULL if none found
inline char* strcut(char* str, const char c)
{
  if (!str) return NULL;
  char* ptr=str;
  while (*ptr!='\0' && *ptr!=c) ptr++;
  if (*ptr=='\0') return NULL;
  *ptr='\0';
  ptr++;
  while (*ptr!='\0' && *ptr==c) ptr++;
  return (*ptr=='\0')? NULL:ptr;
}

// Cuts string at first occurence of charcter c, by overwriting it with '\0'.
// Returns pointer to next char not equal c, OR to end of string '\0' if none found
inline char* strcut_(char* str, const char c)
{
  if (!str) return NULL;
  char* ptr=str;
  while (*ptr!='\0' && *ptr!=c) ptr++;
  if (*ptr=='\0') return ptr;
  *ptr='\0';
  ptr++;
  while (*ptr!='\0' && *ptr==c) ptr++;
  return ptr;
}

// Cuts string at first occurence of substr, by overwriting the first letter with '\0'.
// Returns pointer to next char after occurence of substr, OR to NULL if no such char found 
inline char* strcut(char* str, const char* substr)
{
  char* ptr;     //present location in str being compared to substr
  const char* sptr=substr; //present location in substr being compared to substr
  // while not at end of str and not all of substr is matched yet
  while (1)
    {
      for (ptr=str, sptr=substr; *ptr==*sptr && *ptr!='\0';  ptr++, sptr++) ;
      if (*sptr=='\0') {*str='\0'; return ptr;}
      if (*ptr=='\0')  return NULL;
      str++;
    }
}

// Cuts string at first occurence of substr, by overwriting the first letter with '\0'.
// Returns pointer to next char after occurence of substr, OR to end of string '\0' if no such char found 
inline char* strcut_(char* str, const char* substr)
{
  char* ptr;         //present location in str being compared to substr
  const char* sptr=substr; //present location in substr being compared to str
  // while not at end of str and not all of substr is matched yet
  while (1)
    {
      for (ptr=str, sptr=substr; *ptr==*sptr && *ptr!='\0';  ptr++, sptr++) ;
      if (*sptr=='\0') {*str='\0'; return ptr;}
      if (*ptr=='\0')  return ptr;
      str++;
    }
}

// Copies first word in ptr to str. In other words, copies first block of non whitespace characters, 
// beginning at ptr, to str. If a word is found, returns address of second word in ptr or, if no second
// word is found, returns address to end of word ('\0' character) in ptr string. If no word is found 
// in ptr NULL is returned.
inline char* strwrd(char* str, char* ptr)
{
  ptr=strscn(ptr);    // advance to beginning of next word
  if (ptr) 
    {
      while (*ptr!='\0' && *ptr>32) *(str++) = *(ptr++);
      *str='\0';
      while (*ptr!='\0' && *ptr<=32) ptr++;
      return ptr;
    }
  else return NULL;
}

// Copies first word ***delimited by char c*** in ptr to str. In other words, copies first block of non-c characters, 
// beginning at ptr, to str. If a word is found, returns address of second word in ptr or, if no second
// word is found, returns address to end of word ('\0' character) in ptr string. If no word is found 
// in ptr NULL is returned.
inline char* strwrd(char* str, char* ptr, const char c)
{
  ptr=strscn(ptr,c);    // advance to beginning of next word
  if (ptr) 
    {
      while (*ptr!='\0' && *ptr!=c) *(str++) = *(ptr++);
      *str='\0';
      while (*ptr!='\0' && *ptr==c) ptr++;
     return ptr;
    }
  else return NULL;
}

// Similar to Perl's tr/abc/ABC/: Replaces all chars in str found in one list with characters from the second list
// Returns the number of replaced charactrs
int strtr(char* str, const char oldchars[], const char newchars[])
{
  char* ptr;
  const char *plist;
  int ntr=0;
  for (ptr=str; *ptr!='\0'; ptr++)
    for (plist=oldchars; *plist!='\0'; plist++)
      if (*ptr==*plist) 
	{
	  *ptr=newchars[plist-oldchars]; 
	  ntr++;
	  break;
	}
  return ntr;
}

// Similar to Perl's tr/abc//d: deletes all chars in str found in the list
// Returns number of removed characters
int strtrd(char* str, const char chars[])
{
  char* ptr0=str;
  char* ptr1=str;
  const char *plist;
  while (*ptr1!='\0')
    {
      for (plist=chars; *plist!='\0'; plist++)
	if (*ptr1==*plist) break;
      if (*plist=='\0') {*ptr0=*ptr1; ptr0++;}
      ptr1++;
    }
  return ptr1-ptr0;
}

// Similar to Perl's tr/a-z//d: deletes all chars in str found in the list
// Returns number of removed characters
int strtrd(char* str, char char1, char char2)
{
  char* ptr0=str;
  char* ptr1=str;
  while (*ptr1!='\0')
    {
      if (*ptr1>=char1 && *ptr1<=char2) {*ptr0=*ptr1; ptr0++;}
      ptr1++;
    }
  return ptr1-ptr0;
}

// transforms str into an all uppercase string
char* uprstr(char* str)
{
  char* s=str;
  while (*s !='\0') {if (*s>='a' && *s<='z') *s+='A'-'a';s++;}
  return(str);
}

// transforms str into an all uppercase string
char* lwrstr(char* str)
{
  char* s=str;
  while (*s !='\0') {if (*s>='A' && *s<='Z') *s+='a'-'A'; s++;}
  return(str);
}

// transforms chr into an uppercase character
inline char uprchr(char chr)
{
  return (chr>='a' && chr<='z')? chr+'A'-'a' : chr;
}

// transforms chr into an lowercase character
inline char lwrchr(char chr)
{
  return (chr>='A' && chr<='Z')? chr-'A'+'a' : chr;
}


// Replaces first occurence of str1 by str2 in str. Returns pointer to first occurence or NULL if not found 
// ATTENTION: if str2 is longer than str1, allocated memory of str must be long enough!!
inline char* strsubst(char* str, const char str1[], const char str2[])
{
  char* ptr = strstr(str,str1); 
  strcpy(ptr,str2);
  return ptr;
}

// Gives elapsed time since first call to this function
inline void ElapsedTimeSinceFirstCall(const char str[]) 
{
  timeval t;
  static double tfirst=0;
  if (tfirst==0) 
    {
      gettimeofday(&t, NULL);
      tfirst = 1E-6*t.tv_usec + t.tv_sec;
    }
  gettimeofday(&t, NULL);
  printf("Elapsed time since first call:%12.3fs %s\n",1E-6*t.tv_usec + t.tv_sec - tfirst,str);
}

// Gives elapsed time since last call to this function
inline void ElapsedTimeSinceLastCall(const char str[]) 
{
  timeval t;
  static double tlast=0.0;
  if (tlast==0.0) 
    {
      gettimeofday(&t, NULL);
      tlast = 1.0E-6*t.tv_usec + t.tv_sec;
    }
  gettimeofday(&t, NULL);
  printf("Elapsed time since last call:%12.3fs %s\n",1.0E-6*t.tv_usec + t.tv_sec - tlast,str);
  tlast = 1.0E-6*t.tv_usec + t.tv_sec;
}

inline char* RemovePath(char outname[], char filename[])
{
  char* ptr;
#ifdef WINDOWS
  ptr=strrchr(filename,92);  //return adress for LAST \ (backslash) in name
#else
  ptr=strrchr(filename,'/'); //return adress for LAST / in name
#endif
  if (!ptr) ptr=filename; else ptr++;
  strcpy(outname,ptr);
  return outname;
}

inline char* RemoveExtension(char outname[], char filename[])
{
  char *ptr1;
  ptr1=strrchr(filename,'.');       //return adress for LAST '.' in name
  if (ptr1) {*ptr1='\0'; strcpy(outname,filename); *ptr1='.';} else strcpy(outname,filename);  
  return outname;
}

inline char* RemovePathAndExtension(char outname[], char filename[])
{
  char *ptr, *ptr1;
#ifdef WINDOWS
  ptr=strrchr(filename,92);  //return adress for LAST \ (backslash) in name
#else
  ptr=strrchr(filename,'/'); //return adress for LAST / in name
#endif
  if (!ptr) ptr=filename; else ptr++;
  ptr1=strrchr(filename,'.');       //return adress for LAST '.' in name
  if (ptr1) {*ptr1='\0'; strcpy(outname,ptr); *ptr1='.';} else strcpy(outname,ptr);  
  return outname;
}

inline char* Extension(char extension[], char filename[])
{
  char* ptr;
  ptr=strrchr(filename,'.');      //return adress for LAST '.' in name
  if (ptr) strcpy(extension,ptr+1); else *extension='\0';
  return extension;
}

// Path includes last '/'
inline char* Pathname(char pathname[], char filename[])
{
  char* ptr;
  char chr;
#ifdef WINDOWS
  ptr=strrchr(filename,92);  //return adress for LAST \ (backslash) in name
#else
  ptr=strrchr(filename,'/'); //return adress for LAST / in name
#endif
  if (ptr) {chr=*(++ptr); *ptr='\0'; strcpy(pathname,filename); *ptr=chr;} else *pathname='\0';
  return pathname;
}

// Swaps two integer elements in array k
inline void swapi(int k[], int i, int j)
{
  int temp;
  temp=k[i]; k[i]=k[j]; k[j]=temp;
}

// QSort sorting routine. time complexity of O(N ln(N)) on average
// Sorts the index array k between elements i='left' and i='right' in such a way that afterwards 
// v[k[i]] is sorted downwards (up=-1) or upwards (up=+1)
void QSortInt(int v[], int k[], int left, int right, int up=+1)
{
  int i;      
  int last;   // last element to have been swapped
  
  if (left>=right) return;        // do nothing if less then 2 elements to sort
  // Put pivot element in the middle of the sort range to the side (to position 'left') ...
  swapi(k,left,(left+right)/2);  
  last=left; 
  // ... and swap all elements i SMALLER than the pivot 
  // with an element that is LARGER than the pivot (element last+1):
  if (up==1)
    {
      for (i=left+1; i<=right; i++)
	if (v[k[i]]<v[k[left]]) swapi(k,++last,i);
    }
  else
    for (i=left+1; i<=right; i++)
      if (v[k[i]]>v[k[left]]) swapi(k,++last,i);

  // Put the pivot to the right of the elements which are SMALLER, left to elements which are LARGER
  swapi(k,left,last);

  // Sort the elements left from the pivot and right from the pivot
  QSortInt(v,k,left,last-1,up);
  QSortInt(v,k,last+1,right,up);
}

// QSort sorting routine. time complexity of O(N ln(N)) on average
// Sorts the index array k between elements i='left' and i='right' in such a way that afterwards 
// v[k[i]] is sorted downwards (up=-1) or upwards (up=+1)
void QSortFloat(float v[], int k[], int left, int right, int up=+1)
{
  int i;      
  int last;   // last element to have been swapped
  void swapi(int k[], int i, int j);
  
  if (left>=right) return;        // do nothing if less then 2 elements to sort
  // Put pivot element in the middle of the sort range to the side (to position 'left') ...
  swapi(k,left,(left+right)/2);  
  last=left; 
  // ... and swap all elements i SMALLER than the pivot 
  // with an element that is LARGER than the pivot (element last+1):
  if (up==1)
    {
    for (i=left+1; i<=right; i++)
      if (v[k[i]]<v[k[left]]) swapi(k,++last,i);
    }
  else
    for (i=left+1; i<=right; i++)
      if (v[k[i]]>v[k[left]]) swapi(k,++last,i);

  // Put the pivot to the right of the elements which are SMALLER, left to elements which are LARGER
  swapi(k,left,last);

  // Sort the elements left from the pivot and right from the pivot
  QSortFloat(v,k,left,last-1,up);
  QSortFloat(v,k,last+1,right,up);
}

/**
 * @brief comparison function for qsort,
 * sorts floating point numbers ascendingly
 *
 * @param cv1 ponter to 1st entry to be sorted
 * @param cv2 ponter to 2nd entry to be sorted
 *
 * @return 0 if entries are equal,
 *  +/-1 if 1st greater/smaller than 2nd
 */
int CompFltAsc(const void *cv1, const void *cv2){

  float f1 = *(float *)cv1;
  float f2 = *(float *)cv2;

  if      (f1 > f2) { return +1; }
  else if (f1 < f2) { return -1; }
  else              { return  0; }

} /* this is the end of CompFltAsc() */

//Return random number in the range [0,1]
inline float frand() { return rand()/(RAND_MAX+1.0); }


/////////////////////////////////////////////////////////////////////////////////////
//// Execute system command
/////////////////////////////////////////////////////////////////////////////////////
void runSystem(std::string cmd, int v = 2)
{
  if (v>2)
    cout << "Command: " << cmd << "!\n";
  int res = system(cmd.c_str());
  if (res!=0) 
    {
      cerr << endl << "ERROR when executing: " << cmd << "!\n";
      throw 1;
    }
    
}
