%{

/**************************************************************************
 * Copyright (c) 2019- Guillermo A. Perez
 * 
 * This file is part of the HOA2AIG tool.
 * 
 * HOA2AIG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * HOA2AIG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with HOA2AIG.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Guillermo A. Perez
 * University of Antwerp
 * guillermoalberto.perez@uantwerpen.be
 *************************************************************************/

#include <stdio.h>

#include "hoaparser.h"
%}

%option yylineno

%x C_COMMENT

%%
\/\*                     { BEGIN(C_COMMENT); }
<C_COMMENT>\*\/          { BEGIN(INITIAL); }
<C_COMMENT>.|\n          { }

\"(\.|[^\"])*\"          printf("STRING ");
HOA:                     printf("hoa-header ");
Start:                   printf("start ");
controllable-AP:         printf("controllable-AP ");
AP:                      printf("AP ");
Alias:                   printf("alias ");
Acceptance:              printf("acceptance ");
acc-name:                printf("acc-name ");
tool:                    printf("tool ");
name:                    printf("name ");
properties:              printf("properties ");
State:                   printf("state-header ");
[a-zA-Z_][0-9a-zA-Z_-]*: printf("HEADERNAME ");
0|[1-9][0-9]*            printf("INT ");
[tf]                     printf("BOOLEAN ");
[a-zA-Z_][0-9a-zA-Z_-]*  printf("IDENTIFIER ");
@[0-9a-zA-Z_-]+          printf("ANAME ");
--BODY--                 printf("BEGINBODY ");
--END--                  printf("ENDBODY ");
\(                       printf("LPAR ");
\)                       printf("RPAR ");
\{                       printf("LBRACE ");
\}                       printf("RBRACE ");
\[                       printf("LSQBRACE ");
\]                       printf("RSQBRACE ");
&                        printf("BOOLAND ");
\|                       printf("BOOLOR ");
!                        printf("BOOLNOT ");
[ \t\n\r]+               /* ignore whitespace */;
<<EOF>>                  { printf("(eof %u)\n", yylineno); return 0; }
%%
