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

#include "simplehoa.h"
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
            str, yylineno - 1);  // FIXME: This is shifted for some reason
    autoError = true;
}
%}

/* Yacc declarations: Tokens/terminal used in the grammar */

%locations
%error-verbose

/* HEADER TOKENS */
/* compulsory: must appear exactly once */
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
    NodeType nodetype;
    IntList* numlist;
    StringList* strlist;
    BTree* tree;
}

%token <string> STRING IDENTIFIER ANAME HEADERNAME
%token <number> INT
%token <boolean> BOOL

%type <number> header_item header_list
%type <numlist> state_conj
%type <strlist> string_list id_list boolintid_list
%type <nodetype> accid
%type <tree> acceptance_cond acc_cond_conj acc_cond_atom

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
                  loadedData->version = $2;
                  if (seenHeader[HOAHDR])
                      hdrItemError("HOA:");
                  else
                      seenHeader[HOAHDR] = true;
              };

header_list: /* empty */ { /* no new item, nothing to check */ }
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
           | AP INT string_list                {
                                                 loadedData->noAPs = $2;
                                                 loadedData->aps = $3;
                                                 $$ = AP;
                                               }
           | CNTAP int_list                    { $$ = CNTAP; }
           | ALIAS ANAME label_expr            { $$ = ALIAS; }
           | ACCEPTANCE INT acceptance_cond    { 
                                                 loadedData->noAccSets = $2;
                                                 loadedData->acc $3;
                                                 $$ = ACCEPTANCE;
                                               }
           | ACCNAME IDENTIFIER boolintid_list { 
                                                 loadedData->accNameID = $2;
                                                 loadedData->accNameParameters
                                                    = concatStrLists(
                                                        loadedData->
                                                            accNameParameters,
                                                        $3
                                                    );
                                                 $$ = ACCNAME;
                                               }
           | TOOL STRING maybe_string          { $$ = TOOL; }
           | NAME STRING                       {
                                                 loadedData->name = $2;
                                                 $$ = NAME;
                                               }
           | PROPERTIES id_list                { 
                                                 loadedData->properties =
                                                     concatStrLists(
                                                         loadedData->properties,
                                                         $2
                                                     );
                                                 $$ = PROPERTIES; }
           | HEADERNAME boolintstrid_list      { 
                                                 printf("Headername: %s\n",
                                                        $1);
                                                 $$ = HEADERNAME;
                                               }
           ;

state_conj: INT                 { $$ = newIntNode($1); }
          | state_conj "&" INT  { $$ = appendIntNode($1, $3); }
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

acceptance_cond: acc_cond_conj                     { $$ = $1; }
               | acceptance_cond "|" acc_cond_conj { $$ = orBTree($1, $3); }
               ;
               
acc_cond_conj: acc_cond_atom                   { $$ = $1; }
             | acc_cond_conj "&" acc_cond_atom { $$ = andBTree($1, $3); }
             ;
             
acc_cond_atom: accid "(" INT ")"       { $$ = idBTree($1, $3, false); }
             | accid "(" "!" INT ")"   { $$ = idBTree($1, $4, true); }
             | "(" acceptance_cond ")" { $$ = $2; }
             | BOOL                    { $$ = boolBTree($1); }
             ; 

accid: FIN { $$ = NT_FIN; }
     | INF { $$ = NT_INF; }
     ;

boolintid_list: /* empty */               { $$ = NULL; }
              | boolintid_list BOOL       { 
                                            $$ = $2 ? appendStrNode($1, "True")
                                                    : appendStrNode($1, "False");
                                          }
              | boolintid_list INT        {
                                            char buffer[66];
                                            sprintf(buffer, "%d", $2);
                                            $$ = appendStrNode($1, buffer);
                                          }
              | boolintid_list IDENTIFIER { $$ = appendStrNode($1, $2); }
              ;

boolintstrid_list: /* empty */
                 | boolintstrid_list BOOL
                 | boolintstrid_list INT
                 | boolintstrid_list STRING
                 | boolintstrid_list IDENTIFIER;

string_list: /* empty */        { $$ = NULL; }
           | string_list STRING { $$ = appendStrNode($1, $2); }
           ;

id_list: /* empty */        { $$ = NULL; }
       | id_list IDENTIFIER { $$ = appendStrNode($1, $2); }
       ;

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
