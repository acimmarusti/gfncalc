#ifndef READ_INPUTS_H

#define READ_INPUTS_H

void conv_str_num(std::string in , int &num);

void conv_str_num(std::string in , double &num);

void safe_read_num(std::string message , int &num);

void safe_read_num(std::string message , double &num);

void safe_read_char(std::string message , char &resp) ;

void read_input(int &flag , std::string lookfor , std::string lookhere , std::string touse , int &inputnum);

void read_input(int &flag , std::string lookfor , std::string lookhere , std::string touse , double &inputnum);

#endif
