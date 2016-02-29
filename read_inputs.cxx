#include <iostream>
#include <sstream>
#include <string>

using namespace std;

#include "read_inputs.h"

void conv_str_num(string in , int &num)
{
  stringstream conv(in);
  
  conv >> num;
}

void conv_str_num(string in , double &num)
{
  stringstream conv(in);
  
  conv >> num;
}

void safe_read_num(string message , int &num)
{
  string input = "";	
  
  while (true) {
    cout << message;
    getline(cin, input);
    
    //------This code converts from string to number safely------//
    
    stringstream conv(input);
    
    if (conv >> num) break;
    
    cout << "Invalid number, please try again" << endl;
  }
  
  cout << "You entered: " << num << endl << endl;
}

void safe_read_num(string message , double &num)
{
  string input = "";
  
  while (true) {
    cout << message;
    getline(cin, input);
    
    //------This code converts from string to number safely------//
    
    stringstream conv(input);
    
    if (conv >> num) break;
    
    cout << "Invalid number, please try again" << endl;
  }

  cout << "You entered: " << num << endl << endl;
}

void safe_read_char(string message , char & resp)
{
  string input = "" ;	
  
  while (true){
    cout << message ;
    getline(cin, input) ;
    
    if (input.length() == 1) {
      resp = input[0] ;
      
      break ;
    }
    
    cout << "Invalid character, please try again" << endl ;
  }
  
  cout << "You entered: " << resp << endl << endl ;
}

void read_input(int &flag , string lookfor , string lookhere , string touse , double &inputnum)
{
  if (flag == 0)
    
    if (lookfor == lookhere) {
      
      conv_str_num(touse , inputnum);
      
      flag = 1;
    }
}

void read_input(int &flag , string lookfor , string lookhere , string touse , int &inputnum)
{
  if (flag == 0)
    
    if (lookfor == lookhere) {

      conv_str_num(touse , inputnum);
      
      flag = 1;
    }
}
