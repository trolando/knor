/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     HOAHDR = 1,
     ACCEPTANCE = 2,
     STATES = 3,
     AP = 4,
     CNTAP = 5,
     ACCNAME = 6,
     TOOL = 7,
     NAME = 8,
     START = 9,
     ALIAS = 10,
     PROPERTIES = 11,
     LPAR = 258,
     RPAR = 259,
     LBRACE = 260,
     RBRACE = 261,
     LSQBRACE = 262,
     RSQBRACE = 263,
     BOOLOR = 264,
     BOOLAND = 265,
     BOOLNOT = 266,
     STRING = 267,
     INT = 268,
     BOOL = 269,
     IDENTIFIER = 270,
     ANAME = 271,
     HEADERNAME = 272,
     STATEHDR = 273,
     INF = 274,
     FIN = 275,
     BEGINBODY = 276,
     ENDBODY = 277
   };
#endif
/* Tokens.  */
#define HOAHDR 1
#define ACCEPTANCE 2
#define STATES 3
#define AP 4
#define CNTAP 5
#define ACCNAME 6
#define TOOL 7
#define NAME 8
#define START 9
#define ALIAS 10
#define PROPERTIES 11
#define LPAR 258
#define RPAR 259
#define LBRACE 260
#define RBRACE 261
#define LSQBRACE 262
#define RSQBRACE 263
#define BOOLOR 264
#define BOOLAND 265
#define BOOLNOT 266
#define STRING 267
#define INT 268
#define BOOL 269
#define IDENTIFIER 270
#define ANAME 271
#define HEADERNAME 272
#define STATEHDR 273
#define INF 274
#define FIN 275
#define BEGINBODY 276
#define ENDBODY 277




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef int YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

#if ! defined YYLTYPE && ! defined YYLTYPE_IS_DECLARED
typedef struct YYLTYPE
{
  int first_line;
  int first_column;
  int last_line;
  int last_column;
} YYLTYPE;
# define yyltype YYLTYPE /* obsolescent; will be withdrawn */
# define YYLTYPE_IS_DECLARED 1
# define YYLTYPE_IS_TRIVIAL 1
#endif

extern YYLTYPE yylloc;
