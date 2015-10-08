///////////////////////////////////////////////////////////////////////////////
// File: example.cc                                                          //
// Part of the CmdLine library                                               //
//                                                                           //
// Copyright (c) 2007 Gavin Salam                                            //
//                                                                           //
// This program is free software; you can redistribute it and/or modify      //
// it under the terms of the GNU General Public License as published by      //
// the Free Software Foundation; either version 2 of the License, or         //
// (at your option) any later version.                                       //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program; if not, write to the Free Software               //
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA //
//                                                                           //
// $Revision:: 139                                                          $//
// $Date:: 2007-01-23 16:09:23 +0100 (Tue, 23 Jan 2007)                     $//
///////////////////////////////////////////////////////////////////////////////

#include "CmdLine.hh"
#include<iostream>

using namespace std;

int main (int argc, char ** argv) {

  // construct the cmdline object
  CmdLine cmdline(argc, argv);

  // the value<T> template deduces the correct type from the
  // default value for the option (if present)
  double dval = cmdline.value("-d",0.0);
  cout << "dval = " << dval << endl;

  // for options with a default character value, we usually want
  // a string result -- so this must be specified explicitly
  string sval = cmdline.value<string>("-s","not there");
  cout << "sval = " << sval << endl;

  // time stamping is not strictly purpose of cmdline -- but it's
  // useful anyhow
  cout << "Time now = " << cmdline.time_stamp() << endl;

  // get info on machine
  cout << "Command run on " << cmdline.unix_uname() << endl;

  // make sure we've used all options that were provided on command-line
  cmdline.assert_all_options_used();
}
