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
     LPAR = 258,
     RPAR = 259,
     LBRACE = 260,
     RBRACE = 261,
     LSQBRACE = 262,
     RSQBRACE = 263,
     BOOLAND = 264,
     BOOLOR = 265,
     BOOLNOT = 266,
     STRING = 267,
     INT = 268,
     BOOL = 269,
     IDENTIFIER = 270,
     ANAME = 271,
     HEADERNAME = 272,
     HOAHDR = 273,
     STATES = 274,
     AP = 275,
     ALIAS = 276,
     ACCEPTANCE = 277,
     ACCNAME = 278,
     TOOL = 279,
     NAME = 280,
     PROPERTIES = 281,
     STATEHDR = 282,
     INF = 283,
     FIN = 284,
     BEGINBODY = 285,
     ENDBODY = 286,
     CNTAP = 287,
     START = 288
   };
#endif
/* Tokens.  */
#define LPAR 258
#define RPAR 259
#define LBRACE 260
#define RBRACE 261
#define LSQBRACE 262
#define RSQBRACE 263
#define BOOLAND 264
#define BOOLOR 265
#define BOOLNOT 266
#define STRING 267
#define INT 268
#define BOOL 269
#define IDENTIFIER 270
#define ANAME 271
#define HEADERNAME 272
#define HOAHDR 273
#define STATES 274
#define AP 275
#define ALIAS 276
#define ACCEPTANCE 277
#define ACCNAME 278
#define TOOL 279
#define NAME 280
#define PROPERTIES 281
#define STATEHDR 282
#define INF 283
#define FIN 284
#define BEGINBODY 285
#define ENDBODY 286
#define CNTAP 287
#define START 288




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
