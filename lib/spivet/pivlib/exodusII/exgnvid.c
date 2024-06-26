/*
 * Copyright (c) 2005 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
 * retains certain rights in this software.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.  
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
/*****************************************************************************
*
* exgnvid - ex_get_nodal_varid
*
* entry conditions - 
*   input parameters:
*       int     exoid           exodus file id
*
* exit conditions - 
*       int*    varid           array of nodal variable varids
*
* revision history - 
*
*  $Id: exgnvid.c,v 1.5 2007/10/08 15:01:41 gdsjaar Exp $
*
*****************************************************************************/

#include <stdlib.h>
#include "exodusII.h"
#include "exodusII_int.h"

/*
 *  returns the varids for the nodal variables.
 */

int ex_get_nodal_varid(int exoid, int *varid)
{
  int i, dimid, nvarid;
  long num_vars;
  char errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  if ((dimid = ncdimid (exoid, DIM_NUM_NOD_VAR)) == -1) {
    num_vars = 0;
    if (ncerr == NC_EBADDIM)
      return(EX_NOERR);     /* no nodal variables defined */
    else
      {
        exerrval = ncerr;
        sprintf(errmsg,
                "Error: failed to locate nodal variable names in file id %d",
                exoid);
        ex_err("ex_get_nodal_varid",errmsg,exerrval);
        return (EX_FATAL);
      }
  }

  if (ncdiminq (exoid, dimid, (char *) 0, &num_vars) == -1) {
    exerrval = ncerr;
    sprintf(errmsg,
            "Error: failed to get number of nodal variables in file id %d",
            exoid);
    ex_err("ex_get_nodal_varid",errmsg,exerrval);
    return (EX_FATAL);
  }
   
  if (ex_large_model(exoid) == 0) {
    /* All varids are the same; */
    if ((nvarid = ncvarid (exoid, VAR_NOD_VAR)) == -1) {
      exerrval = ncerr;
      sprintf(errmsg,
              "Warning: could not find nodal variables in file id %d",
              exoid);
      ex_err("ex_get_nodal_varid",errmsg,exerrval);
      return (EX_WARN);
    }
    for (i=0; i < num_vars; i++) {
      varid[i] = nvarid;
    }
  } else {
    /* Variables stored separately; each has a unique varid */
    for (i=0; i < num_vars; i++) {
      if ((nvarid = ncvarid (exoid, VAR_NOD_VAR_NEW(i+1))) == -1) {
        exerrval = ncerr;
        sprintf(errmsg,
                "Warning: could not find nodal variable %d in file id %d",
                i+1, exoid);
        ex_err("ex_get_nodal_varid",errmsg,exerrval);
        return (EX_WARN);
      }
      varid[i] = nvarid;
    }
  }
  return(EX_NOERR);
}
