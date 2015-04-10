/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#define MAXARGS 100

typedef struct {
 char *str;
 int *flag;
 int type;
 char **arg;
} CmdLineData;

#define NOARG 0
#define INTARG 1
#define FLTARG 2
#define STRARG 3
#define FILARG 4
#define OPTARG 5




