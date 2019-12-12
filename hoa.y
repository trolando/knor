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
#include <stdbool.h>

#include "hoa.h"
#include "hoalexer.h"

HoaData* loadedData;

void yyerror(const char* str) {
    fprintf(stderr, "Parsing error: %s [line %d]\n", str, yylineno);
}
 
int yywrap() {
    return 1;
}

bool autoError = false;
bool seenHeader[12] = {false};
const char* headerStrs[] = {"", "HOA", "Acceptance", "States", "AP",
                            "controllable-AP", "acc-name", "tool",
                            "name", "Start", "Alias", "properties"};

void hdrItemError(const char* str) {
    fprintf(stderr,
            "Automaton error: more than one %s header item [line %d]\n",
            str, yylineno);
    autoError = true;
}
%}

/* Yacc declarations: Tokens/terminal used in the grammar */

%locations
%error-verbose

/* HEADER TOKENS */
/* compulsory */
%token HOAHDR 1 ACCEPTANCE 2 /* indexed from 1 since 0 means EOF for bison */
/* at most once */
%token STATES 3 AP 4 CNTAP 5 ACCNAME 6 TOOL 7 NAME 8
/* multiple */
%token START 9 ALIAS 10 PROPERTIES 11

/* OTHERS */
%token LPAR "("
%token RPAR ")"
%token LBRACE "{"
%token RBRACE "}"
%token LSQBRACE "["
%token RSQBRACE "]"
%token BOOLOR "|"
%token BOOLAND "&"
%token BOOLNOT "!"
%token STATEHDR INF FIN BEGINBODY ENDBODY

%union {
    int number;
    char* string;
    bool boolean;
    IntList* numlist;
}

%token <string> STRING IDENTIFIER ANAME HEADERNAME
%token <number> INT
%token <boolean> BOOL

%type <number> header_item header_list
%type <numlist> state_conj

%%
/* Grammar rules and actions follow */

automaton: header BEGINBODY body ENDBODY
         {
            if (!seenHeader[HOAHDR]) /* redundant because of the grammar */
                yyerror("No HOA: header item");
            if (!seenHeader[ACCEPTANCE])
                yyerror("No Acceptance: header item");
         };

header: format_version header_list;

format_version: HOAHDR IDENTIFIER
              {
                  if (seenHeader[HOAHDR])
                      hdrItemError("HOA:");
                  else
                      seenHeader[HOAHDR] = true;
              };

header_list: /* empty */
           | header_list header_item
           {
               if ($2 <= 7) {
                   if (seenHeader[$2])
                       hdrItemError(headerStrs[$2]);
                   else
                       seenHeader[$2] = true;
               }
           };

header_item: STATES INT                        { 
                                                 loadedData->noStates = $2;
                                                 $$ = STATES;
                                               }
           | START state_conj                  {
                                                 loadedData->start = $2;
                                                 $$ = START;
                                               }
           | AP INT string_list                { $$ = AP; }
           | CNTAP int_list                    { $$ = CNTAP; }
           | ALIAS ANAME label_expr            { $$ = ALIAS; }
           | ACCEPTANCE INT acceptance_cond    { $$ = ACCEPTANCE; }
           | ACCNAME IDENTIFIER boolintid_list { $$ = ACCNAME; }
           | TOOL STRING maybe_string          { $$ = TOOL; }
           | NAME STRING                       { $$ = NAME; }
           | PROPERTIES id_list                { $$ = PROPERTIES; }
           | HEADERNAME boolintstrid_list      { $$ = HEADERNAME; }
           ;

state_conj: INT                 { $$ = newIntNode($1); }
          | state_conj "&" INT { $$ = appendIntNode($1, $3); }
          ;

label_expr: lab_exp_conj
          | label_expr "|" lab_exp_conj;
          
lab_exp_conj: lab_exp_atom
            | lab_exp_conj "&" lab_exp_atom;

lab_exp_atom: BOOL
            | INT
            | ANAME
            | "!" lab_exp_atom
            | "(" label_expr ")";

acceptance_cond: acc_cond_conj
               | acceptance_cond "|" acc_cond_conj;
               
acc_cond_conj: acc_cond_atom
             | acc_cond_conj "&" acc_cond_atom;
             
acc_cond_atom: accid "(" INT ")"
             | accid "(" "!" INT ")"
             | "(" acceptance_cond ")"
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
  
int parseHoa(FILE* input, HoaData* data) {
    loadedData = data;
    yyin = input;
    int ret = yyparse();
    return ret | autoError;
}
