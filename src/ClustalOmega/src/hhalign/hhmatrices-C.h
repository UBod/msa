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
 * RCS $Id: hhmatrices-C.h 274 2012-04-24 23:28:24Z dave $
 */

// Substitution matrices and their background frequencies

// The following background frequencies were calculated by the formula pa = (P[a,b]/(pa*pb))^(-1) * (1,...,1)
// For the Blousum50-matrix this becomes  pb[a]= SUM_(b=1,20) (2^(BLOSUM50[a,b]/3))^(-1) 
//                     A    R    N    D    C    Q    E    G    H    I    L     K    M    F    P    S    T    W    Y    V 
// Gonnet            7.68,5.14,4.02,5.41,1.89,3.27,5.99,7.56,3.69,5.06,10.01,5.97,2.20,3.50,4.54,4.67,7.12,1.25,3.95,7.28
// BLOSUM50          8.24,6.24,4.46,4.77,2.03,2.90,6.78,6.69,2.53,6.89,10.7 ,5.04,1.49,4.93,3.97,5.95,6.13,1.34,3.45,6.28

const float Gonnet[]={
//  A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V 
 10227, 3430, 2875, 3869, 1625, 2393, 4590, 6500, 2352, 3225, 5819, 4172, 1435, 1579, 3728, 4610, 6264,  418, 1824, 5709, // A
  3430, 7780, 2209, 2589,  584, 2369, 3368, 3080, 2173, 1493, 3093, 5701,  763,  859, 1893, 2287, 3487,  444, 1338, 2356, // R 
  2875, 2209, 3868, 3601,  501, 1541, 2956, 3325, 1951, 1065, 2012, 2879,  532,  688, 1480, 2304, 3204,  219, 1148, 1759, // N
  3869, 2589, 3601, 8618,  488, 2172, 6021, 4176, 2184, 1139, 2151, 3616,  595,  670, 2086, 2828, 3843,  204, 1119, 2015, // D
  1625,  584,  501,  488, 5034,  355,  566,  900,  516,  741, 1336,  591,  337,  549,  419,  901, 1197,  187,  664, 1373, // C
  2393, 2369, 1541, 2172,  355, 1987, 2891, 1959, 1587, 1066, 2260, 2751,  570,  628, 1415, 1595, 2323,  219,  871, 1682, // Q
  4590, 3368, 2956, 6021,  566, 2891, 8201, 3758, 2418, 1624, 3140, 4704,  830,  852, 2418, 2923, 4159,  278, 1268, 2809, // E
  6500, 3080, 3325, 4176,  900, 1959, 3758,26066, 2016, 1354, 2741, 3496,  741,  797, 2369, 3863, 4169,  375, 1186, 2569, // G
  2352, 2173, 1951, 2184,  516, 1587, 2418, 2016, 5409, 1123, 2380, 2524,  600, 1259, 1298, 1642, 2446,  383,  876, 1691, // H
  3225, 1493, 1065, 1139,  741, 1066, 1624, 1354, 1123, 6417, 9630, 1858, 1975, 2225, 1260, 1558, 3131,  417, 1697, 7504, // I
  5819, 3093, 2012, 2151, 1336, 2260, 3140, 2741, 2380, 9630,25113, 3677, 4187, 5540, 2670, 2876, 5272, 1063, 3945,11005, // L
  4172, 5701, 2879, 3616,  591, 2751, 4704, 3496, 2524, 1858, 3677, 7430,  949,  975, 2355, 2847, 4340,  333, 1451, 2932, // K
  1435,  763,  532,  595,  337,  570,  830,  741,  600, 1975, 4187,  949, 1300, 1111,  573,  743, 1361,  218,  828, 2310, // M
  1579,  859,  688,  670,  549,  628,  852,  797, 1259, 2225, 5540,  975, 1111, 6126,  661,  856, 1498, 1000, 4464, 2602, // F
  3728, 1893, 1480, 2086,  419, 1415, 2418, 2369, 1298, 1260, 2670, 2355,  573,  661,11834, 2320, 3300,  179,  876, 2179, // P
  4610, 2287, 2304, 2828,  901, 1595, 2923, 3863, 1642, 1558, 2876, 2847,  743,  856, 2320, 3611, 4686,  272, 1188, 2695, // S
  6264, 3487, 3204, 3843, 1197, 2323, 4159, 4169, 2446, 3131, 5272, 4340, 1361, 1498, 3300, 4686, 8995,  397, 1812, 5172, // T
   418,  444,  219,  204,  187,  219,  278,  375,  383,  417, 1063,  333,  218, 1000,  179,  272,  397, 4101, 1266,  499, // W
  1824, 1338, 1148, 1119,  664,  871, 1268, 1186,  876, 1697, 3945, 1451,  828, 4464,  876, 1188, 1812, 1266, 9380, 2227, // Y
  5709, 2356, 1759, 2015, 1373, 1682, 2809, 2569, 1691, 7504,11005, 2932, 2310, 2602, 2179, 2695, 5172,  499, 2227,11569};// V

const float Blosum30[]= {
//    A      R      N      D      C      Q      E      G      H      I      L      K      M      F      P      S      T      W      Y      V 
  0.0096,
  0.0038,0.0109,
  0.0031,0.0019,0.0055,
  0.0043,0.0026,0.0027,0.0095,
  0.0014,0.0011,0.0010,0.0010,0.0070,
  0.0031,0.0031,0.0014,0.0018,0.0007,0.0039,
  0.0044,0.0031,0.0024,0.0037,0.0018,0.0028,0.0094,
  0.0052,0.0032,0.0032,0.0035,0.0011,0.0020,0.0035,0.0173,
  0.0016,0.0014,0.0011,0.0011,0.0004,0.0010,0.0018,0.0015,0.0060,
  0.0040,0.0022,0.0024,0.0018,0.0012,0.0016,0.0023,0.0036,0.0012,0.0072,
  0.0056,0.0039,0.0032,0.0043,0.0023,0.0027,0.0051,0.0051,0.0022,0.0066,0.0139,
  0.0044,0.0043,0.0027,0.0032,0.0010,0.0021,0.0053,0.0039,0.0013,0.0026,0.0044,0.0063,
  0.0018,0.0012,0.0009,0.0008,0.0004,0.0007,0.0012,0.0013,0.0008,0.0014,0.0027,0.0017,0.0012,
  0.0027,0.0023,0.0017,0.0011,0.0008,0.0010,0.0016,0.0022,0.0008,0.0027,0.0055,0.0023,0.0008,0.0077,
  0.0028,0.0021,0.0012,0.0021,0.0007,0.0015,0.0030,0.0030,0.0014,0.0017,0.0027,0.0029,0.0005,0.0011,0.0091,
  0.0056,0.0035,0.0028,0.0034,0.0013,0.0021,0.0038,0.0051,0.0017,0.0033,0.0047,0.0042,0.0010,0.0027,0.0024,0.0075,
  0.0040,0.0019,0.0024,0.0024,0.0009,0.0016,0.0023,0.0028,0.0010,0.0027,0.0046,0.0025,0.0010,0.0016,0.0020,0.0041,0.0046,
  0.0005,0.0007,0.0002,0.0004,0.0003,0.0004,0.0007,0.0011,0.0002,0.0005,0.0009,0.0006,0.0002,0.0007,0.0004,0.0005,0.0003,0.0027,
  0.0014,0.0022,0.0008,0.0015,0.0004,0.0011,0.0016,0.0016,0.0010,0.0018,0.0045,0.0017,0.0007,0.0024,0.0011,0.0017,0.0014,0.0009,0.0044,
  0.0056,0.0031,0.0022,0.0027,0.0012,0.0015,0.0026,0.0033,0.0012,0.0063,0.0074,0.0030,0.0015,0.0032,0.0017,0.0036,0.0036,0.0005,0.0024,0.0083};

const float Blosum40[]= {
//    A      R      N      D      C      Q      E      G      H      I      L      K      M      F      P      S      T      W      Y      V 
  0.0148,
  0.0029,0.0109,
  0.0029,0.0019,0.0069,
  0.0032,0.0021,0.0031,0.0126,
  0.0015,0.0007,0.0007,0.0009,0.0093,
  0.0026,0.0024,0.0017,0.0016,0.0004,0.0048,
  0.0040,0.0026,0.0023,0.0048,0.0010,0.0032,0.0118,
  0.0066,0.0023,0.0035,0.0031,0.0012,0.0020,0.0028,0.0260,
  0.0014,0.0012,0.0013,0.0014,0.0003,0.0010,0.0015,0.0014,0.0060,
  0.0037,0.0017,0.0017,0.0016,0.0008,0.0012,0.0018,0.0024,0.0009,0.0105,
  0.0050,0.0030,0.0021,0.0027,0.0015,0.0023,0.0036,0.0034,0.0017,0.0082,0.0209,
  0.0041,0.0051,0.0027,0.0030,0.0010,0.0026,0.0045,0.0035,0.0013,0.0023,0.0037,0.0099,
  0.0017,0.0010,0.0007,0.0007,0.0003,0.0007,0.0010,0.0013,0.0007,0.0018,0.0037,0.0012,0.0018,
  0.0022,0.0016,0.0013,0.0013,0.0008,0.0008,0.0015,0.0021,0.0009,0.0031,0.0055,0.0018,0.0011,0.0105,
  0.0027,0.0014,0.0014,0.0019,0.0005,0.0013,0.0026,0.0028,0.0008,0.0020,0.0021,0.0023,0.0007,0.0010,0.0151,
  0.0060,0.0025,0.0030,0.0030,0.0013,0.0026,0.0034,0.0052,0.0013,0.0025,0.0034,0.0034,0.0011,0.0019,0.0022,0.0089,
  0.0040,0.0019,0.0022,0.0024,0.0011,0.0015,0.0025,0.0029,0.0009,0.0028,0.0039,0.0029,0.0011,0.0019,0.0022,0.0043,0.0070,
  0.0007,0.0005,0.0003,0.0003,0.0001,0.0004,0.0005,0.0008,0.0002,0.0005,0.0009,0.0005,0.0002,0.0008,0.0003,0.0004,0.0003,0.0045,
  0.0019,0.0016,0.0010,0.0011,0.0004,0.0010,0.0014,0.0017,0.0012,0.0020,0.0031,0.0017,0.0010,0.0032,0.0009,0.0015,0.0014,0.0008,0.0060,
  0.0054,0.0023,0.0017,0.0022,0.0012,0.0014,0.0025,0.0028,0.0008,0.0083,0.0082,0.0027,0.0018,0.0031,0.0018,0.0032,0.0040,0.0005,0.0020,0.0113};
 
const float Blosum50[]= {
//    A      R      N      D      C      Q      E      G      H      I      L      K      M      F      P      S      T      W      Y      V 
  0.0192,
  0.0027,0.0152,
  0.0024,0.0020,0.0101,
  0.0026,0.0019,0.0035,0.0161,
  0.0015,0.0005,0.0006,0.0005,0.0091,
  0.0022,0.0025,0.0016,0.0017,0.0004,0.0057,
  0.0034,0.0029,0.0023,0.0048,0.0006,0.0033,0.0141,
  0.0062,0.0020,0.0031,0.0028,0.0009,0.0017,0.0023,0.0316,
  0.0012,0.0013,0.0015,0.0011,0.0003,0.0010,0.0013,0.0011,0.0064,
  0.0035,0.0015,0.0013,0.0012,0.0008,0.0011,0.0015,0.0018,0.0007,0.0140,
  0.0048,0.0028,0.0017,0.0018,0.0014,0.0019,0.0026,0.0027,0.0013,0.0104,0.0304,
  0.0033,0.0064,0.0027,0.0026,0.0006,0.0029,0.0043,0.0028,0.0014,0.0017,0.0027,0.0130,
  0.0016,0.0009,0.0007,0.0005,0.0004,0.0008,0.0008,0.0009,0.0005,0.0022,0.0042,0.0010,0.0029,
  0.0020,0.0012,0.0009,0.0008,0.0006,0.0007,0.0012,0.0015,0.0009,0.0030,0.0058,0.0012,0.0012,0.0154,
  0.0022,0.0011,0.0011,0.0015,0.0004,0.0011,0.0019,0.0019,0.0006,0.0013,0.0017,0.0018,0.0005,0.0007,0.0171,
  0.0062,0.0025,0.0032,0.0028,0.0011,0.0022,0.0030,0.0044,0.0012,0.0021,0.0029,0.0031,0.0010,0.0016,0.0018,0.0111,
  0.0039,0.0021,0.0026,0.0022,0.0010,0.0015,0.0024,0.0025,0.0009,0.0029,0.0038,0.0026,0.0011,0.0015,0.0016,0.0047,0.0100,
  0.0005,0.0004,0.0002,0.0002,0.0001,0.0003,0.0004,0.0005,0.0002,0.0005,0.0008,0.0004,0.0003,0.0009,0.0002,0.0003,0.0004,0.0059,
  0.0015,0.0013,0.0009,0.0009,0.0004,0.0009,0.0012,0.0012,0.0013,0.0018,0.0027,0.0013,0.0007,0.0039,0.0006,0.0013,0.0012,0.0008,0.0077,
  0.0054,0.0020,0.0015,0.0016,0.0013,0.0014,0.0020,0.0022,0.0007,0.0107,0.0092,0.0022,0.0021,0.0030,0.0016,0.0029,0.0041,0.0005,0.0018,0.0164};

const float Blosum65[]= {
//  A      R      N      D      C      Q      E      G      H      I      L      K      M      F      P      S      T      W      Y      V
  0.0222,
  0.0022,0.0181,
  0.0019,0.0019,0.0148,
  0.0021,0.0015,0.0037,0.0225,
  0.0016,0.0004,0.0004,0.0004,0.0127,
  0.0018,0.0024,0.0015,0.0016,0.0003,0.0076,
  0.0029,0.0025,0.0021,0.0049,0.0003,0.0034,0.0168,
  0.0057,0.0016,0.0027,0.0024,0.0007,0.0013,0.0018,0.0396,
  0.0010,0.0013,0.0014,0.0009,0.0002,0.0010,0.0013,0.0009,0.0096,
  0.0031,0.0012,0.0009,0.0011,0.0012,0.0008,0.0012,0.0013,0.0006,0.0191,
  0.0043,0.0023,0.0013,0.0014,0.0016,0.0016,0.0019,0.0020,0.0009,0.0115,0.0388,
  0.0032,0.0062,0.0024,0.0024,0.0005,0.0030,0.0040,0.0024,0.0012,0.0015,0.0024,0.0166,
  0.0013,0.0007,0.0005,0.0004,0.0004,0.0007,0.0006,0.0007,0.0003,0.0025,0.0052,0.0008,0.0045,
  0.0016,0.0009,0.0007,0.0007,0.0005,0.0005,0.0008,0.0011,0.0008,0.0030,0.0056,0.0009,0.0012,0.0186,
  0.0021,0.0009,0.0008,0.0011,0.0003,0.0008,0.0014,0.0013,0.0005,0.0009,0.0013,0.0015,0.0004,0.0005,0.0195,
  0.0065,0.0022,0.0030,0.0027,0.0011,0.0018,0.0029,0.0037,0.0011,0.0017,0.0024,0.0030,0.0008,0.0012,0.0016,0.0137,
  0.0037,0.0017,0.0022,0.0019,0.0009,0.0013,0.0021,0.0022,0.0007,0.0027,0.0033,0.0023,0.0010,0.0012,0.0013,0.0048,0.0133,
  0.0004,0.0003,0.0002,0.0001,0.0002,0.0002,0.0002,0.0004,0.0002,0.0004,0.0007,0.0003,0.0002,0.0009,0.0001,0.0003,0.0003,0.0074,
  0.0013,0.0009,0.0007,0.0006,0.0004,0.0006,0.0008,0.0008,0.0016,0.0015,0.0023,0.0010,0.0006,0.0043,0.0004,0.0010,0.0009,0.0010,0.0113,
  0.0049,0.0015,0.0011,0.0012,0.0014,0.0011,0.0016,0.0017,0.0006,0.0120,0.0094,0.0019,0.0023,0.0025,0.0012,0.0023,0.0035,0.0004,0.0015,0.0206};


const float Blosum80[]= {
//  A      R      N      D      C      Q      E      G      H      I      L      K      M      F      P      S      T      W      Y      V
  0.0252,
  0.0020,0.0210,
  0.0016,0.0017,0.0166,
  0.0018,0.0013,0.0037,0.0262,
  0.0015,0.0003,0.0004,0.0003,0.0172,
  0.0017,0.0024,0.0014,0.0014,0.0003,0.0094,
  0.0028,0.0023,0.0019,0.0048,0.0003,0.0035,0.0208,
  0.0053,0.0015,0.0025,0.0022,0.0006,0.0011,0.0017,0.0463,
  0.0009,0.0012,0.0012,0.0008,0.0002,0.0011,0.0012,0.0008,0.0104,
  0.0027,0.0010,0.0007,0.0008,0.0011,0.0007,0.0010,0.0009,0.0004,0.0220,
  0.0036,0.0018,0.0011,0.0011,0.0014,0.0014,0.0015,0.0016,0.0008,0.0111,0.0442,
  0.0029,0.0061,0.0022,0.0020,0.0004,0.0028,0.0036,0.0020,0.0010,0.0012,0.0019,0.0190,
  0.0011,0.0006,0.0004,0.0003,0.0004,0.0007,0.0006,0.0005,0.0003,0.0025,0.0052,0.0007,0.0053,
  0.0014,0.0007,0.0006,0.0006,0.0005,0.0005,0.0006,0.0009,0.0007,0.0027,0.0052,0.0007,0.0010,0.0211,
  0.0021,0.0009,0.0007,0.0009,0.0003,0.0007,0.0012,0.0010,0.0004,0.0007,0.0012,0.0012,0.0003,0.0004,0.0221,
  0.0064,0.0020,0.0029,0.0024,0.0010,0.0017,0.0026,0.0034,0.0010,0.0015,0.0021,0.0026,0.0007,0.0010,0.0014,0.0167,
  0.0036,0.0015,0.0020,0.0016,0.0009,0.0012,0.0019,0.0019,0.0007,0.0024,0.0028,0.0020,0.0009,0.0011,0.0011,0.0048,0.0156,
  0.0003,0.0002,0.0001,0.0001,0.0001,0.0002,0.0002,0.0003,0.0001,0.0003,0.0006,0.0002,0.0002,0.0007,0.0001,0.0002,0.0002,0.0087,
  0.0011,0.0007,0.0006,0.0005,0.0003,0.0005,0.0006,0.0006,0.0016,0.0013,0.0020,0.0008,0.0005,0.0046,0.0003,0.0009,0.0008,0.0010,0.0148,
  0.0046,0.0013,0.0009,0.0010,0.0013,0.0010,0.0015,0.0014,0.0005,0.0123,0.0089,0.0015,0.0022,0.0022,0.0010,0.0021,0.0033,0.0004,0.0012,0.0246};

const float mism = 2000; // standard mismatch score
const float matc = 7000;

// RNA matrix, Gonnet type format  
const float ribosum[]={
//  A     C     G     T     U     R     Y     M     K     S     W     B     D     H     V     N     ?     ?     ?     ? 
  9500,  600,  700,  900,  900, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // A
   600, 4000,  300, 1200, 1200, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // C 
   700,  300, 3500,  700,  700, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // G
   900, 1200,  700, 8200, 8200, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // T
   900, 1200,  700, 8200, 8200, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // U
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // R
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // Y
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // M
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // K
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // S
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // W
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // B
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // D
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // H
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // V
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // N
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // ?
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // ?
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // ?
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism};// ?

const float dna_basic[]={
//  A     C     G     T     U     R     Y     M     K     S     W     B     D     H     V     N     ?     ?     ?     ? 
  matc, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // A
  mism, matc, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // C 
  mism, mism, matc, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // G
  mism, mism, mism, matc, matc, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // T
  mism, mism, mism, matc, matc, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // U
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // R
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // Y
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // M
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // K
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // S
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // W
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // B
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // D
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // H
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // V
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // N
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // ?
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // ?
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // ?
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism};// ?
  
const float dave_nucleo[]={
//  A     C     G     T     U     R     Y     M     K     S     W     B     D     H     V     N     ?     ?     ?     ? 
  8000, 1500, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // A
  1500, 6500, 1000, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // C 
  mism, 1000, 6500, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // G
  mism, mism, mism, 7000, 7000, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // T
  mism, mism, mism, 7000, 7000, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // U
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // R
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // Y
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // M
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // K
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // S
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // W
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // B
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // D
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // H
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // V
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // N
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // ?
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // ?
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, // ?
  mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism, mism};// ?
  
// prediction accuracy of Psipred: 
// Ppred[cf][B][A] = P(A,B,cf)/P(A)/P(B,cf) = P(A|B,cf)/P(A)
// A = observed ss state  B = predicted ss state  cf = confidence value of prediction 
//float Ppred[MAXCF][NSSPRED][NDSSP]= 
const float Ppred[]= 
  {
//pred/obs  -      H      E      ~      S      T      G      B
        1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=-
        1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // H
        1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // E
        1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // ~
        1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=0
        1.000, 1.128, 0.519, 0.834, 0.957, 1.488, 2.106, 1.085 , // H
        1.000, 0.233, 2.240, 1.216, 0.913, 0.519, 0.923, 1.759 , // E
        1.000, 0.640, 1.017, 1.122, 1.069, 1.242, 2.140, 1.999 , // ~
        1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=1
        1.000, 1.251, 0.485, 0.771, 0.847, 1.371, 2.266, 0.864 , // H
        1.000, 0.222, 2.542, 1.069, 0.804, 0.428, 0.671, 1.728 , // E
        1.000, 0.474, 1.103, 1.295, 1.232, 1.214, 1.835, 1.989 , // ~
        1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=2
        1.000, 1.383, 0.426, 0.637, 0.778, 1.349, 2.436, 0.824 , // H
        1.000, 0.202, 2.769, 0.999, 0.714, 0.320, 0.551, 1.566 , // E
        1.000, 0.395, 1.005, 1.407, 1.376, 1.336, 1.725, 2.063 , // ~
        1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=3
        1.000, 1.531, 0.369, 0.552, 0.682, 1.280, 2.420, 0.698 , // H
        1.000, 0.169, 2.970, 0.954, 0.556, 0.273, 0.489, 1.504 , // E
        1.000, 0.352, 0.843, 1.515, 1.542, 1.456, 1.684, 1.958 , // ~
        1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=4
        1.000, 1.750, 0.305, 0.444, 0.537, 1.134, 2.295, 0.600 , // H
        1.000, 0.124, 3.179, 0.847, 0.513, 0.228, 0.413, 1.897 , // E
        1.000, 0.282, 0.718, 1.664, 1.630, 1.577, 1.625, 1.877 , // ~
        1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=5
        1.000, 1.952, 0.250, 0.353, 0.456, 0.982, 2.050, 0.466 , // H
        1.000, 0.102, 3.464, 0.699, 0.453, 0.174, 0.284, 1.357 , // E
        1.000, 0.227, 0.574, 1.782, 1.846, 1.681, 1.418, 1.885 , // ~
        1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=6
        1.000, 2.183, 0.171, 0.263, 0.319, 0.792, 1.933, 0.345 , // H
        1.000, 0.079, 3.712, 0.612, 0.281, 0.133, 0.196, 1.089 , // E
        1.000, 0.173, 0.458, 1.915, 2.007, 1.766, 1.220, 1.704 , // ~
        1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=7
        1.000, 2.389, 0.132, 0.192, 0.224, 0.605, 1.605, 0.183 , // H
        1.000, 0.053, 3.997, 0.449, 0.201, 0.072, 0.141, 0.919 , // E
        1.000, 0.109, 0.328, 2.013, 2.304, 1.882, 0.960, 1.512 , // ~
        1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=8
        1.000, 2.668, 0.065, 0.098, 0.144, 0.354, 1.059, 0.102 , // H
        1.000, 0.029, 4.285, 0.284, 0.113, 0.044, 0.059, 0.522 , // E
        1.000, 0.053, 0.200, 2.099, 2.444, 2.133, 0.671, 1.290 , // ~
        1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=9
        1.000, 2.966, 0.009, 0.023, 0.036, 0.113, 0.214, 0.017 , // H
        1.000, 0.010, 4.555, 0.119, 0.027, 0.010, 0.013, 0.209 , // E
        1.000, 0.026, 0.101, 2.576, 1.853, 2.204, 0.308, 0.499   // ~
  };

// float Ppred[]= 
//   {
// //pred/obs  -      H      E      ~      S      T      G      B
//         1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=-
//         1.000, 2.000, 0.500, 1.000, 1.000, 1.000, 1.000, 0.500,  // H
//         1.000, 0.500, 2.000, 0.500, 1.000, 1.000, 0.500, 1.000,  // E
//         1.000, 1.000, 1.000, 2.000, 2.000, 2.000, 1.000, 1.000,  // ~
//         1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=0
//         1.000, 2.000, 0.500, 1.000, 1.000, 1.000, 1.000, 0.500,  // H
//         1.000, 0.500, 2.000, 0.500, 1.000, 1.000, 0.500, 1.000,  // E
//         1.000, 1.000, 1.000, 2.000, 2.000, 2.000, 1.000, 1.000,  // ~
//         1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=1
//         1.000, 2.000, 0.500, 1.000, 1.000, 1.000, 1.000, 0.500,  // H
//         1.000, 0.500, 2.000, 0.500, 1.000, 1.000, 0.500, 1.000,  // E
//         1.000, 1.000, 1.000, 2.000, 2.000, 2.000, 1.000, 1.000,  // ~
//         1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=2
//         1.000, 2.000, 0.500, 1.000, 1.000, 1.000, 1.000, 0.500,  // H
//         1.000, 0.500, 2.000, 0.500, 1.000, 1.000, 0.500, 1.000,  // E
//         1.000, 1.000, 1.000, 2.000, 2.000, 2.000, 1.000, 1.000,  // ~
//         1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=3
//         1.000, 2.000, 0.500, 1.000, 1.000, 1.000, 1.000, 0.500,  // H
//         1.000, 0.500, 2.000, 0.500, 1.000, 1.000, 0.500, 1.000,  // E
//         1.000, 1.000, 1.000, 2.000, 2.000, 2.000, 1.000, 1.000,  // ~
//         1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=4
//         1.000, 2.000, 0.500, 1.000, 1.000, 1.000, 1.000, 0.500,  // H
//         1.000, 0.500, 2.000, 0.500, 1.000, 1.000, 0.500, 1.000,  // E
//         1.000, 1.000, 1.000, 2.000, 2.000, 2.000, 1.000, 1.000,  // ~
//         1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=5
//         1.000, 2.000, 0.500, 1.000, 1.000, 1.000, 1.000, 0.500,  // H
//         1.000, 0.500, 2.000, 0.500, 1.000, 1.000, 0.500, 1.000,  // E
//         1.000, 1.000, 1.000, 2.000, 2.000, 2.000, 1.000, 1.000,  // ~
//         1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=6
//         1.000, 2.000, 0.500, 1.000, 1.000, 1.000, 1.000, 0.500,  // H
//         1.000, 0.500, 2.000, 0.500, 1.000, 1.000, 0.500, 1.000,  // E
//         1.000, 1.000, 1.000, 2.000, 2.000, 2.000, 1.000, 1.000,  // ~
//         1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=7
//         1.000, 2.000, 0.500, 1.000, 1.000, 1.000, 1.000, 0.500,  // H
//         1.000, 0.500, 2.000, 0.500, 1.000, 1.000, 0.500, 1.000,  // E
//         1.000, 1.000, 1.000, 2.000, 2.000, 2.000, 1.000, 1.000,  // ~
//         1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=8
//         1.000, 2.000, 0.500, 1.000, 1.000, 1.000, 1.000, 0.500,  // H
//         1.000, 0.500, 2.000, 0.500, 1.000, 1.000, 0.500, 1.000,  // E
//         1.000, 1.000, 1.000, 2.000, 2.000, 2.000, 1.000, 1.000,  // ~
//         1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 , // - conf=9
//         1.000, 2.000, 0.500, 1.000, 1.000, 1.000, 1.000, 0.500,  // H
//         1.000, 0.500, 2.000, 0.500, 1.000, 1.000, 0.500, 1.000,  // E
//         1.000, 1.000, 1.000, 2.000, 2.000, 2.000, 1.000, 1.000,  // ~


float Pobs[]={0, 0.3268,0.2119,0.2061,0.0913,0.1143,0.0376,0.0120};

float Pevo_full[]=
  {
//     -     H     E     C    
    1.00, 0.00, 0.00, 0.00,
    0.00, 0.94, 0.00, 0.04,  
    0.00, 0.00, 0.92, 0.04,
    0.00, 0.06, 0.08, 0.92
  };

//psipred accuracy for confidence values 0-9
const float p_acc[]={0.00,0.47,0.53,0.56,0.58,0.62,0.69,0.74,0.82,0.88,0.96}; 

/**
 * @brief
 */
void 
SetBlosumMatrix(const float BlosumXX[])
{
  int a,b,n=0;
  if (v>=3) printf("Using the BLOSUM%2i matrix\n",par.matrix);
  for (a=0; a<20; ++a)
    for (pb[a]=0.0f, b=0; b<=a; ++b,++n)
      P[a][b] = BlosumXX[n];
  for (a=0; a<19; a++)
    for (b=a+1; b<20; ++b)
      P[a][b] = P[b][a];
  for (a=0; a<20; ++a) P[a][20]=P[20][a]=1.0f;
  return;
}

/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Set (global variable) substitution matrix with derived matrices and background frequencies
 */
void 
SetSubstitutionMatrix()
{
  int a,b;
  switch (par.matrix)
    {
    default:
    case 0:  //Gonnet matrix
      if (v>=3) cout<<"Using the Gonnet matrix ";
      for (a=0; a<20; ++a)
	for (pb[a]=0.0f, b=0; b<20; ++b)
	  P[a][b] = 0.000001f*Gonnet[a*20+b];
      for (a=0; a<20; ++a) P[a][20]=P[20][a]=1.0f;
      break;

    case 30:  //BLOSUM30
      SetBlosumMatrix(Blosum30);
      break;
    case 40:  //BLOSUM40
      SetBlosumMatrix(Blosum40);
      break;
    case 50:  //BLOSUM50
      SetBlosumMatrix(Blosum50);
      break;
    case 65:  //BLOSUM65
      SetBlosumMatrix(Blosum65);
      break;
    case 80:  //BLOSUM80
      SetBlosumMatrix(Blosum80);
      break;
   }
  
  // Check transition probability matrix, renormalize P and calculate pb[a]
  float sumab=0.0f;
  for (a=0; a<20; a++)
    for (b=0; b<20; ++b) sumab+=P[a][b];
  for (a=0; a<20; a++)
    for (b=0; b<20; ++b) P[a][b]/=sumab;
  for (a=0; a<20; a++)
    for (pb[a]=0.0f, b=0; b<20; ++b) pb[a]+=P[a][b];

  //Compute similarity matrix for amino acid pairs (for calculating consensus sequence)
  for (a=0; a<20; ++a)
    for (b=0; b<20; ++b)
      Sim[a][b] = P[a][b]*P[a][b]/P[a][a]/P[b][b];

  //Precompute matrix R for amino acid pseudocounts:
  for (a=0; a<20; ++a)
    for (b=0; b<20; ++b)   
      R[a][b] = P[a][b]/pb[b]; //R[a][b]=P(a|b)
  
  //Precompute matrix R for amino acid pseudocounts:
  for (a=0; a<20; ++a)
    for (b=0; b<20; ++b)   
      S[a][b] = log2(R[a][b]/pb[a]); // S[a][b] = log2(P(a,b)/P(a)/P(b))
  
  // Evaluate sequence identity underlying substitution matrix
  if (v>=3)
    {
      float id=0.0f;
      float entropy=0.0f; 
      float entropy_pb=0.0f;
      float mut_info=0.0f;
      for (a=0; a<20; ++a) id+=P[a][a];
      for (a=0; a<20; ++a) entropy_pb-=pb[a]*log2(pb[a]);
      for (a=0; a<20; ++a) 
	  for (b=0; b<20; ++b) 
	    {
	      entropy-=P[a][b]*log2(R[a][b]);
	      mut_info += P[a][b]*S[a][b];
	    }
      
      printf(": sequence identity = %2.0f%%; entropy per column = %4.2f bits (out of %4.2f); mutual information = %4.2f bits\n",100*id,entropy,entropy_pb,mut_info);
    }

  if (v>=4) //Debugging: probability matrix and dissimilarity matrix 
    {
      cout<<"Check matrix: before renormalization sum P(a,b)= "<<sumab<<"...\n";//PRINT
      cout<<"      A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
      cout<<"p[] ";
      for (a=0; a<20; a++)  printf("%4.1f ",100*pb[a]);
      cout<<endl<<"\nSubstitution matrix log2( P(a,b)/p(a)/p(b) ) (in bits):\n";
      cout<<"      A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
      for (b=0; b<20; b++)
	{
	  cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%4.1f ",S[a][b]);
	  cout<<endl;
	}
      cout<<endl<<"\nOdds matrix P(a,b)/p(a)/p(b):\n";
      cout<<"      A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
      for (b=0; b<20; b++)
	{
	  cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%4.1f ",P[b][a]/pb[a]/pb[b]);
	  cout<<endl;
	}
      cout<<endl<<"\nMatrix of conditional probabilities P(a|b) = P(a,b)/p(b) (in %):\n";
      cout<<"      A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
      for (b=0; b<20; b++)
	{
	  cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%4.1f ",100*R[b][a]);
	  cout<<endl;
	}
      cout<<endl<<"\nProbability matrix P(a,b) (in %):\n";
      cout<<"      A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V\n";
      for (b=0; b<20; b++)
	{
	  cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%5.0f ",1000000*P[b][a]);
	  cout<<endl;
	}
      cout<<endl<<"Similarity matrix P(a,b)^2/P(a,a,)/P(b,b) (in %):\n";
      cout<<"      A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
      for (b=0; b<20; b++)
	{
	  cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%4.0f ",100*Sim[b][a]);
	  cout<<endl;
	}
      cout<<endl;


    }
}
 
/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Set (global variable) substitution matrix with derived matrices and background frequencies
 */
void 
SetRnaSubstitutionMatrix()
{
	printf("SET RNA SUBSTITUTION MATRIX ....");
  int a,b;
  for (a=0; a<20; ++a)
	for (pb[a]=0.0f, b=0; b<20; ++b)
	  P[a][b] = 0.000001f*ribosum[a*20+b];
      for (a=0; a<20; ++a) P[a][20]=P[20][a]=1.0f;
  
  // Check transition probability matrix, renormalize P and calculate pb[a]
  float sumab=0.0f;
  for (a=0; a<20; a++)
    for (b=0; b<20; ++b) sumab+=P[a][b];
  for (a=0; a<20; a++)
    for (b=0; b<20; ++b) P[a][b]/=sumab;
  for (a=0; a<20; a++)
    for (pb[a]=0.0f, b=0; b<20; ++b) pb[a]+=P[a][b];

  //Compute similarity matrix for amino acid pairs (for calculating consensus sequence)
  for (a=0; a<20; ++a)
    for (b=0; b<20; ++b)
      Sim[a][b] = P[a][b]*P[a][b]/P[a][a]/P[b][b];

  //Precompute matrix R for amino acid pseudocounts:
  for (a=0; a<20; ++a)
    for (b=0; b<20; ++b)   
      R[a][b] = P[a][b]/pb[b]; //R[a][b]=P(a|b)
  
  //Precompute matrix R for amino acid pseudocounts:
  for (a=0; a<20; ++a)
    for (b=0; b<20; ++b)   
      S[a][b] = log2(R[a][b]/pb[a]); // S[a][b] = log2(P(a,b)/P(a)/P(b))
  
  // Evaluate sequence identity underlying substitution matrix
  if (v>=3)
    {
      float id=0.0f;
      float entropy=0.0f; 
      float entropy_pb=0.0f;
      float mut_info=0.0f;
      for (a=0; a<20; ++a) id+=P[a][a];
      for (a=0; a<20; ++a) entropy_pb-=pb[a]*log2(pb[a]);
      for (a=0; a<20; ++a) 
	  for (b=0; b<20; ++b) 
	    {
	      entropy-=P[a][b]*log2(R[a][b]);
	      mut_info += P[a][b]*S[a][b];
	    }
      
      printf(": sequence identity = %2.0f%%; entropy per column = %4.2f bits (out of %4.2f); mutual information = %4.2f bits\n",100*id,entropy,entropy_pb,mut_info);
    }

  if (v>=4) //Debugging: probability matrix and dissimilarity matrix 
    {
      cout<<"Check matrix: before renormalization sum P(a,b)= "<<sumab<<"...\n";//PRINT
      cout<<"      A    C    G    T    U    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?\n";
      cout<<"p[] ";
      for (a=0; a<20; a++)  printf("%4.1f ",100*pb[a]);
      cout<<endl<<"\nSubstitution matrix log2( P(a,b)/p(a)/p(b) ) (in bits):\n";
      cout<<"      A    C    G    T    U    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?\n";
      for (b=0; b<20; b++)
	{
	  cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%4.1f ",S[a][b]);
	  cout<<endl;
	}
      cout<<endl<<"\nOdds matrix P(a,b)/p(a)/p(b):\n";
      cout<<"      A    C    G    T    U    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?\n";
      for (b=0; b<20; b++)
	{
	  cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%4.1f ",P[b][a]/pb[a]/pb[b]);
	  cout<<endl;
	}
      cout<<endl<<"\nMatrix of conditional probabilities P(a|b) = P(a,b)/p(b) (in %):\n";
      cout<<"      A    C    G    T    U    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?\n";
      for (b=0; b<20; b++)
	{
	  cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%4.1f ",100*R[b][a]);
	  cout<<endl;
	}
      cout<<endl<<"\nProbability matrix P(a,b) (in %):\n";
      cout<<"      A    C    G    T    U    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?\n";
      for (b=0; b<20; b++)
	{
	  cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%5.0f ",1000000*P[b][a]);
	  cout<<endl;
	}
      cout<<endl<<"Similarity matrix P(a,b)^2/P(a,a,)/P(b,b) (in %):\n";
      cout<<"      A    C    G    T    U    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?\n";
      for (b=0; b<20; b++)
	{
	  cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%4.0f ",100*Sim[b][a]);
	  cout<<endl;
	}
      cout<<endl;


    }
}
 
/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Set (global variable) substitution matrix with derived matrices and background frequencies
 */
void 
SetDnaSubstitutionMatrix()
{

  int a,b;
  for (a=0; a<20; ++a)
	for (pb[a]=0.0f, b=0; b<20; ++b)
	  P[a][b] = 0.000001f*dna_basic[a*20+b];
      for (a=0; a<20; ++a) P[a][20]=P[20][a]=1.0f;
  
  // Check transition probability matrix, renormalize P and calculate pb[a]
  float sumab=0.0f;
  for (a=0; a<20; a++)
    for (b=0; b<20; ++b) sumab+=P[a][b];
  for (a=0; a<20; a++)
    for (b=0; b<20; ++b) P[a][b]/=sumab;
  for (a=0; a<20; a++)
    for (pb[a]=0.0f, b=0; b<20; ++b) pb[a]+=P[a][b];

  //Compute similarity matrix for amino acid pairs (for calculating consensus sequence)
  for (a=0; a<20; ++a)
    for (b=0; b<20; ++b)
      Sim[a][b] = P[a][b]*P[a][b]/P[a][a]/P[b][b];

  //Precompute matrix R for amino acid pseudocounts:
  for (a=0; a<20; ++a)
    for (b=0; b<20; ++b)   
      R[a][b] = P[a][b]/pb[b]; //R[a][b]=P(a|b)
  
  //Precompute matrix R for amino acid pseudocounts:
  for (a=0; a<20; ++a)
    for (b=0; b<20; ++b)   
      S[a][b] = log2(R[a][b]/pb[a]); // S[a][b] = log2(P(a,b)/P(a)/P(b))
  
  // Evaluate sequence identity underlying substitution matrix
  if (v>=3)
    {
      float id=0.0f;
      float entropy=0.0f; 
      float entropy_pb=0.0f;
      float mut_info=0.0f;
      for (a=0; a<20; ++a) id+=P[a][a];
      for (a=0; a<20; ++a) entropy_pb-=pb[a]*log2(pb[a]);
      for (a=0; a<20; ++a) 
	  for (b=0; b<20; ++b) 
	    {
	      entropy-=P[a][b]*log2(R[a][b]);
	      mut_info += P[a][b]*S[a][b];
	    }
      
      printf(": sequence identity = %2.0f%%; entropy per column = %4.2f bits (out of %4.2f); mutual information = %4.2f bits\n",100*id,entropy,entropy_pb,mut_info);
    }

  if (v>=4) //Debugging: probability matrix and dissimilarity matrix 
    {
      cout<<"Check matrix: before renormalization sum P(a,b)= "<<sumab<<"...\n";//PRINT
      cout<<"      A    C    G    T    U    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?\n";
      cout<<"p[] ";
      for (a=0; a<20; a++)  printf("%4.1f ",100*pb[a]);
      cout<<endl<<"\nSubstitution matrix log2( P(a,b)/p(a)/p(b) ) (in bits):\n";
      cout<<"      A    C    G    T    U    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?\n";
      for (b=0; b<20; b++)
	{
	  cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%4.1f ",S[a][b]);
	  cout<<endl;
	}
      cout<<endl<<"\nOdds matrix P(a,b)/p(a)/p(b):\n";
      cout<<"      A    C    G    T    U    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?\n";
      for (b=0; b<20; b++)
	{
	  cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%4.1f ",P[b][a]/pb[a]/pb[b]);
	  cout<<endl;
	}
      cout<<endl<<"\nMatrix of conditional probabilities P(a|b) = P(a,b)/p(b) (in %):\n";
      cout<<"      A    C    G    T    U    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?\n";
      for (b=0; b<20; b++)
	{
	  cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%4.1f ",100*R[b][a]);
	  cout<<endl;
	}
      cout<<endl<<"\nProbability matrix P(a,b) (in %):\n";
      cout<<"      A    C    G    T    U    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?\n";
      for (b=0; b<20; b++)
	{
	  cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%5.0f ",1000000*P[b][a]);
	  cout<<endl;
	}
      cout<<endl<<"Similarity matrix P(a,b)^2/P(a,a,)/P(b,b) (in %):\n";
      cout<<"      A    C    G    T    U    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?    ?\n";
      for (b=0; b<20; b++)
	{
	  cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%4.0f ",100*Sim[b][a]);
	  cout<<endl;
	}
      cout<<endl;


    }
}

/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Set secondary structure substitution matrix
 */
void 
SetSecStrucSubstitutionMatrix()
{
  int A;        //observed ss state (determined dssp)
  int B,BB;     //predicted ss states (by psipred)
  int cf,ccf;   //confidence value of prediction
  float P73[NDSSP][NSSPRED][MAXCF];  //P73[cf][B][A] = P(A,B,cf)/P(A)/P(B,cf) = P(A|B,cf)/P(A)
  float sum;

  // S73[A][B][cf][b] = score for matching observed ss state A in query with state B in template
  // predicted with confidence cf, when query and template columns are diverged by b units 
  for (cf=0; cf<MAXCF; cf++)
    for (A=0; A<NDSSP; A++)
      for (B=0; B<NSSPRED; B++)
	{
	  P73[A][B][cf] = 1.-par.ssa + par.ssa*Ppred[cf*NSSPRED*NDSSP + B*NDSSP + A];
	  S73[A][B][cf] = log2(P73[A][B][cf]);
	}

  for (B=0; B<NSSPRED; B++)
    for (cf=0; cf<MAXCF; cf++)
      for (BB=0; BB<NSSPRED; BB++)
	for (ccf=0; ccf<MAXCF; ccf++)
	  {
	    sum=0;
	    for (A=1; A<NDSSP; A++)
	      sum += P73[A][B][cf] * P73[A][BB][ccf] * Pobs[A];
	    S33[B][cf][BB][ccf] = log2(sum);
	  }  
} /* this is the end of SetSecStrucSubstitutionMatrix() */



/*
 * EOF hhmatrices-C.h
 */
