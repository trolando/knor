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

/* C declarations */

#include <stdio.h>
#include <string.h>

#include "hoalexer.h"

void yyerror(const char *str) {
    fprintf(stderr, "Parsing error: %s [line %d]\n", str, yylineno);
}
 
int yywrap() {
    return 1;
} 
%}

/* Yacc declarations: Tokens/terminal used in the grammar */

%locations
%error-verbose

%token LPAR "("
%token RPAR ")"
%token LBRACE "{"
%token RBRACE "}"
%token LSQBRACE "["
%token RSQBRACE "]"
%token BOOLAND "&"
%token BOOLOR "|"
%token BOOLNOT "!"
%token STRING INT BOOL IDENTIFIER ANAME HEADERNAME HOAHDR
%token STATES AP ALIAS ACCEPTANCE ACCNAME TOOL NAME PROPERTIES
%token STATEHDR INF FIN BEGINBODY ENDBODY CNTAP START

%%
/* Grammar rules and actions follow */

automaton: header BEGINBODY body ENDBODY;

header: format_version header_list;

format_version: HOAHDR IDENTIFIER;

header_list: /* empty */
           | header_list header_item;

header_item: STATES INT
           | START state_conj
           | AP INT string_list
           | CNTAP int_list
           | ALIAS ANAME label_expr
           | ACCEPTANCE INT acceptance_cond
           | ACCNAME IDENTIFIER boolintid_list
           | TOOL STRING maybe_string
           | NAME STRING
           | PROPERTIES id_list
           | HEADERNAME boolintstrid_list;

state_conj: INT
          | state_conj "&" INT;

label_expr: BOOL
          | INT
          | ANAME
          | "!" label_expr
          | "(" label_expr ")"
          | label_expr "&" label_expr
          | label_expr "|" label_expr;

acceptance_cond: accid "(" INT ")"
               | accid "(" "!" INT ")"
               | "(" acceptance_cond ")"
               | acceptance_cond "&" acceptance_cond
               | acceptance_cond "|" acceptance_cond
               | BOOL;

accid: FIN
     | INF;

boolintid_list: /* empty */
              | boolintid_list BOOL
              | boolintid_list INT
              | boolintid_list IDENTIFIER;

boolintstrid_list: /* empty */
                 | boolintstrid_list BOOL
                 | boolintstrid_list INT
                 | boolintstrid_list STRING
                 | boolintstrid_list IDENTIFIER;

string_list: /* empty */
           | string_list STRING;

id_list: /* empty */
       | id_list IDENTIFIER;

body: statespec_list;

statespec_list: /* empty */
              | statespec_list state_name trans_list;

state_name: STATEHDR maybe_label INT maybe_string maybe_accsig;

maybe_label: /* empty */
           | "[" label_expr "]";

maybe_string: /* empty */
            | STRING;

maybe_accsig: /* empty */
            | "{" int_list "}";

int_list: /* empty */
        | int_list INT;

trans_list: /* empty */
          | trans_list maybe_label state_conj maybe_accsig;

%%
/* Additional C code */
  
int main() {
    yyparse();
}
