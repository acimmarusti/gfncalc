#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <vector>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#define chunk 20000
#define apd_time_res 164.61e-12

using namespace std;

//---Other header files---//

#include "read_inputs.h"

void press_enter_continue();


int main() {
  cout << endl << "Correlation function calculator V1.2.0" << endl;
  cout << "Written by A. Cimmarusti, D. Norris, E. Cahoon and B. Patterson"
       << endl << endl;

  //---Configuration file and input variable declaration & initialization---//

  int file_counter, first_file_number = 0, last_file_number = 0;
  int bin, window, config_exist;
  int ch_apda, ch_apdb, ch_apdc, ch_off_pulse;
  double pgen_delay, edge_delay, tol_delay;

  string directory = "";
  string input = "";
  string useful_in = "";

  char answer = {0};
  char option = {0};


  //---Reading inputs from configuration file---//

  fstream config_file;

  config_file.open("gfncalc.conf" , ios::in);

  if (config_file.is_open()) {

    while (getline(config_file , input)) {

      size_t ini = input.find_first_not_of(" \t\n\v");

      size_t eq_pos = input.find_last_of("=");

      if (ini != string::npos && input[ini] == '#') continue;

      useful_in = input.substr(eq_pos + 1);

      size_t data = input.find("data");

      size_t apda = input.find("apda");

      size_t apdb = input.find("apdb");

      size_t apdc = input.find("apdc");

      size_t binsize = input.find("binsize");

      size_t win = input.find("window");

      size_t off_pulse = input.find("pulse");

      size_t pgendelay = input.find("pgen");

      size_t toldelay = input.find("tolerance");

      size_t edgedelay = input.find("edge");

      if (data != string::npos) directory = useful_in;

      if (apda != string::npos) conv_str_num(useful_in , ch_apda);

      if (apdb != string::npos) conv_str_num(useful_in , ch_apdb);

      if (apdc != string::npos) conv_str_num(useful_in , ch_apdc);

      if (binsize != string::npos) conv_str_num(useful_in , bin);

      if (win != string::npos) conv_str_num(useful_in , window);

      if (off_pulse != string::npos) conv_str_num(useful_in , ch_off_pulse);

      if (pgendelay != string::npos) conv_str_num(useful_in , pgen_delay);

      if (toldelay != string::npos) conv_str_num(useful_in , tol_delay);

      if (edgedelay != string::npos) conv_str_num(useful_in , edge_delay);

    }

    config_exist = 1;

  } else {

    config_exist = 0;

    cout << "WARNING: Configuration file not present, or could not be read"
	 << endl << "the program will create one" << endl;
  }

  config_file.close();


  //---Asking for inputs---//

  do {

    if (config_exist == 1) {

      cout << "Configuration file settings are:" << endl;

      cout << "d - data path : " << directory << endl;

      cout << "a - channel APD a : " << ch_apda << endl;

      cout << "b - channel APD b : " << ch_apdb << endl;

      cout << "c - channel APD c : " << ch_apdc << endl;

      cout << "w - bin size of width : " << bin
	   << " ( " << bin * apd_time_res << " s ) " << endl;

      cout << "t - window size : " << window
	   << " ( " << window * apd_time_res << " s ) " << endl;

      cout << "p - channel off pulse : " << ch_off_pulse << endl;

      cout << "i - Internal (cable + generator) delay : " << pgen_delay
	   << " ( " << pgen_delay * apd_time_res << " s ) " << endl;

      cout << "e - Delay tolerance : " << tol_delay 
	   << " ( " << tol_delay * apd_time_res << " s ) " << endl;

      cout << "r - Rising edge delay (s) : " << edge_delay << endl;

      safe_read_char("Press the appropriate character to change the desired input or g to continue\n> " , answer);
    }

    if (answer == 'd' || config_exist == 0) {

      directory = "";

      cout << "type in the address of the directory where the data is held (e.g. C:\\data\\ or /home/user/data/)\n> ";

      getline(cin , directory);

      cout << "You entered: " << directory << endl << endl;
    }

    if (answer == 'a' || config_exist == 0)
      safe_read_num("channel APD a : " , ch_apda);

    if (answer == 'b' || config_exist == 0)
      safe_read_num("channel APD b : " , ch_apdb);

    if (answer == 'c' || config_exist == 0)
      safe_read_num("channel APD c : " , ch_apdc);

    if (answer == 'w' || config_exist == 0)
      safe_read_num("bin size of width : " , bin);

    if (answer == 't' || config_exist == 0)
      safe_read_num("window size : " , window);

    if (answer == 'p' || config_exist == 0)
      safe_read_num("channel off pulse : " , ch_off_pulse);

    if (answer == 'i' || config_exist == 0)
      safe_read_num("Internal (cable + generator) delay : " , pgen_delay);

    if (answer == 'e' || config_exist == 0)
      safe_read_num("Delay tolerance : " , tol_delay);

    if (answer == 'r' || config_exist == 0)
      safe_read_num("Rising edge delay (s) : " , edge_delay);

    config_exist = 1;
  }
  while (answer != 'g');


  //---Analysis selection---//

  do {

    safe_read_char("Select analysis to perform:\na - g2\nb - pulsed g2\n> " ,
		   option);
    //    safe_read_char("Select analysis to perform:\na - g2\nb - pulsed g2\nc - g2 and g3\nd - g3\n> " , option);

  }
  while (option != 'a' && option != 'b');
  //  while (option != 'a' && option != 'b' && option != 'c' && option != 'd' );

  string type_analysis;

  if (option == 'a')

    type_analysis = "_g2";

  else if (option == 'b')

    type_analysis = "_g2";


  //---Data selection---//

  cout << "WARNING: the data must be named as data#.asc" << endl;

  safe_read_num("begin with what data file?: " , first_file_number);

  safe_read_num("end with what data file?: " , last_file_number);


  //---Creation/updating of configuration file---//

  config_file.open("gfncalc.conf" , ios::out | ios::trunc);

  config_file << "# Configuration file for histogram calculator\n";

  config_file << "datapath="<< directory << "\n";

  config_file << "apda_ch="<< ch_apda << "\n";

  config_file << "apdb_ch="<< ch_apdb << "\n";

  config_file << "apdc_ch="<< ch_apdc << "\n";

  config_file << "binsize="<< bin << "\n";

  config_file << "window="<< window << "\n";

  config_file << "off_pulse_ch="<< ch_off_pulse << "\n";

  config_file << "pgen_delay="<< pgen_delay << "\n";

  config_file << "tolerance="<< tol_delay << "\n";

  config_file << "edge_delay="<< edge_delay << "\n";

  config_file.close();

  cout << "this program will sequentially analyze data"
       << first_file_number << " to data" << last_file_number << endl;

  press_enter_continue();


  //---Total delay estimation---//

  edge_delay = edge_delay / apd_time_res;

  double total_delay = pgen_delay + edge_delay;


  //---File Operations and iteration---//

  for (file_counter = first_file_number;
       file_counter <= last_file_number;
       file_counter ++) {

    stringstream filename , filenamefinaloutput;

    filename << directory << "data" << file_counter << ".asc";

    filenamefinaloutput << directory << "data" << file_counter
			<< type_analysis << ".asc";

    cout << "opening file " << filename.str() << endl;

    ifstream ifp;

    ifp.open(filename.str().c_str());

    if (!ifp.is_open()) {
      cout << endl << "the file " << filename.str() << " could not be opened"
	   << endl;

      continue;
    }


    //---File variable declaration & initialization---//

    time_t start, end;

    int runcount, stopper, range, array_sub, numbins, chunknum = 0;
    int total_counts_apda = 0, total_counts_apdb = 0, total_counts_apdc = 0;
    int total_counts_pulse = 0, total_triggers = 0;

    int total_valid_clicks = 0;

    double first_time_apda, first_time_apdb, first_time_apdc, first_time_pulse;
    double diff, total_time, click_time, true_total_time, average_n, timing;
    double rate_apda, rate_apdb, rate_apdc, rate_pulse;


    //---Histogram creation and initialization---//

    // Full range [- window , window] plus zero!
    numbins = (int) 2 * window / bin + 1;

    double * g2_times = new (nothrow) double [numbins];

    int * g2_counts = new (nothrow) int [numbins];

    if (g2_times == 0 || g2_counts == 0) {

      cout << "Error: memory could not be allocated to store the g2 counts"
	   << endl;

      exit(1);
    }

    //---g2 histogram---//

    for (runcount = 0; runcount < numbins; runcount ++) {

      g2_times[runcount] = ((runcount * bin) - window);

      g2_counts[runcount] = 0;
    }

    cout << "g2 histogram and error array have been created." << endl;

    /*
    //---g2 and g3 histograms---//

    for (runcount = 0; runcount < numbins; runcount ++)
    {
    g2_times[runcount] = ((runcount * bin) - window);

    g2_counts[runcount] = 0;

    for (g3counter = 0; g3counter < numbins; g3counter ++)
    {
    g3_times[runcount * numbins + g3counter] = ((runcount * numbins + g3counter) * bin - window);

    g3_counts[runcount * numbins + g3counter] = 0;
    }
    }

    cout << "g2 and g3 histogram and error array has been created." << endl;

    //---g3 histogram---//

    for (runcount = 0; runcount < numbins; runcount ++)

    for (g3counter = 0; g3counter < numbins; g3counter ++)
    {
    g3_times[runcount * numbins + g3counter] = ((runcount * numbins + g3counter) * bin - window);

    g3_counts[runcount * numbins + g3counter] = 0;
    }

    cout << "g2 and g3 histogram and error array has been created." << endl;
    */

    time(& start);

    //---Channel/APD sorting---//

    while (!ifp.eof()) {

      //---Chunk variable declaration & initialization---//

      int chunksize, click_channel;
      int chunk_apda_size = 0, chunk_apdb_size = 0, chunk_apdc_size = 0;
      int chunk_pulse_size = 0, chunk_triggers = 0;

      int chunk_valid_clicks = 0;

      double chunk_array_apda[chunk], chunk_array_apdb[chunk];
      double chunk_array_apdc[chunk], chunk_array_pulse[chunk];

      int valid_click_marker[chunk];


      //---Sorting in chunks of predefined number of clicks---//
      for (chunksize = 0; chunksize < chunk; chunksize ++) {

	ifp >> click_time;

	ifp >> click_channel;

	if (ifp.eof()) break;

	if (click_channel == ch_apda) {

	  chunk_array_apda[chunk_apda_size] = click_time;

	  if (chunk_apda_size == 0 && chunknum == 0)
	    first_time_apda = click_time;

	  chunk_apda_size ++;
	}

	if (click_channel == ch_apdb) {

	  chunk_array_apdb[chunk_apdb_size] = click_time;

	  if (chunk_apdb_size == 0 && chunknum == 0)
	    first_time_apdb = click_time;

	  chunk_apdb_size ++;
	}

	if (click_channel == ch_apdc) {

	  chunk_array_apdc[chunk_apdc_size] = click_time;

	  if (chunk_apdc_size == 0 && chunknum == 0)
	    first_time_apdc = click_time;

	  chunk_apdc_size ++;
	}

	if (click_channel == ch_off_pulse) {

	  chunk_array_pulse[chunk_pulse_size] = click_time;

	  if (chunk_pulse_size == 0 && chunknum == 0)
	    first_time_pulse = click_time;

	  chunk_pulse_size ++;
	}

	valid_click_marker[chunksize] = 1;

      }

      //---Turn-off pulses---//

      //---Array containers declaration & initialization---//

      if (option == 'b') {

	/*APD clicks that trigger a turn-off pulse*/

	double pulse_lag;

	int slider = 0;

	for (int pulsepos = 0;  pulsepos < chunk_pulse_size; pulsepos++) {

	  for (int apd_pos = slider; apd_pos < chunk_apda_size; apd_pos++) {

	    pulse_lag = chunk_array_pulse[pulsepos] - chunk_array_apda[apd_pos];

	    if (pulse_lag > total_delay - tol_delay &&
		pulse_lag < total_delay + tol_delay) {

	      valid_click_marker[apd_pos] = 1;

	      chunk_triggers++;

	      slider = apd_pos + 1;

	      break;

	    } else
	      valid_click_marker[apd_pos] = 0;
	  }

	}

	for (int ii = 0; ii < slider; ii++)

	  chunk_valid_clicks += valid_click_marker[ii];

      }

      total_counts_apda += chunk_apda_size;

      total_counts_apdb += chunk_apdb_size;

      total_counts_apdc += chunk_apdc_size;

      total_counts_pulse += chunk_pulse_size;

      total_triggers += chunk_triggers;

      total_valid_clicks += chunk_valid_clicks;

      chunknum ++;


      //---Filling g2 histogram---//

      /* The idea is to find the clicks in apdb that are within 2*window of
	 the clicks in apda. However, since we want the times associated with
	 the histogram bins to be centered. For example center bin with range
	 (- bin/2 , bin/2] , thus having time zero centered This requires two
	 important changes :
	 1. The time associated with index array_sub is shifted by - bin/2
	 The result of this, is that ranges become (] and not [)
	 The latter case is perfectly valid but requires the floor function
	 and the shift would be bin/2. Also adjustments to the comparison 
	 operators when checking clicks
	 2. The window is expanded by bin/2 on both limits
	 This solves the problem of having only half a bin associated with 
	 the boundaries of the array */

      stopper = 0;

      runcount = 0;

      while (runcount < chunk_apda_size) {

	if (valid_click_marker[runcount] == 1) {

	  while (stopper < chunk_apdb_size) {

	    if (chunk_array_apdb[stopper] >
		chunk_array_apda[runcount] - window - bin / 2)
	      break;

	    stopper ++;
	  }

	  range = stopper;

	  while (range < chunk_apdb_size &&
		 chunk_array_apdb[range] <=
		 chunk_array_apda[runcount] + window + bin / 2) {

	    diff = chunk_array_apdb[range] - chunk_array_apda[runcount];

	    array_sub = ceil ((double)(diff + window) / bin - 0.5);

	    //This precaution protects against erroneous memory allocations
	    if (array_sub >= 0 && array_sub < numbins)

	      g2_counts[array_sub]++;

	    else cout << "WARNING : counts out of range" << endl;

	    range ++;
	  }

	}
	runcount++;
      }

      cout << "counting completed " << chunksize << "-click-size chunk "
	   << chunknum << " in data" << file_counter << endl;
    }

    time(& end);

    ifp.close();

    cout << "the histogram array has been completed" << endl;


    //---Number of clicks---//

    cout << "APD a - " << total_counts_apda << " clicks." << endl;

    cout << "APD b - " << total_counts_apdb << " clicks." << endl;

    cout << "APD c - " << total_counts_apdc << " clicks." << endl;

    if (option == 'b') {

      cout << "Off pulses - " << total_counts_pulse << " clicks." << endl;

      cout << "Triggers - " << total_triggers << " clicks." << endl;

      cout << "Valid clicks - " << total_valid_clicks << " clicks." << endl;

    }

    //---Total time data length---//

    if (first_time_apda < first_time_apdb)
      total_time = click_time - first_time_apda;

    else
      total_time = click_time - first_time_apdb;

    true_total_time = total_time * apd_time_res;

    cout.precision(4);

    cout << "the data spans " << fixed << true_total_time << " seconds" << endl;

    timing = difftime(end , start);

    cout << "Processing time for data" << file_counter << " " << timing
	 << " seconds" << endl << endl;

    cout << "calculating g2's and errors" << endl;


    //---Rates and average---//

    if (option == 'b')

      average_n = (double) total_triggers * total_counts_apdb * bin / total_time;
    else

      average_n = (double) total_counts_apda * total_counts_apdb * bin / total_time;

    rate_apda = total_counts_apda / true_total_time;

    rate_apdb = total_counts_apdb / true_total_time;

    rate_apdc = total_counts_apdc / true_total_time;

    rate_pulse = total_counts_pulse / true_total_time;


    //---Output file operations---//

    ofstream finaloutputfp;

    finaloutputfp.open(filenamefinaloutput.str().c_str() , ios::trunc);

    finaloutputfp << "# data " << file_counter << " final output\n";

    finaloutputfp << "# total clicks on APD A\t" << total_counts_apda << "\t\t";

    finaloutputfp << "rate on APD A\t" << rate_apda << "\n";

    finaloutputfp << "# total clicks on APD B\t" << total_counts_apdb << "\t\t";

    finaloutputfp << "rate on APD B\t" << rate_apdb << "\n";

    finaloutputfp << "# total clicks on APD C\t" << total_counts_apdc << "\t\t";

    finaloutputfp << "rate on APD C\t" << rate_apdc << "\n";

    finaloutputfp << "# total number of pulses\t" << total_counts_pulse
		  << "\t\t";

    finaloutputfp << "rate on pulses\t" << rate_pulse << "\n";

    finaloutputfp << "# <n>\t" << average_n << "\n\n";

    finaloutputfp << "# t(us)\t";

    finaloutputfp << "count\t";

    finaloutputfp << "g2\t";

    finaloutputfp << "error\n";


    //---Histogram with g2 and errors---//

    for (runcount = 0; runcount < numbins; runcount ++) {
      finaloutputfp << setprecision(6)
		    << g2_times[runcount] * apd_time_res * 1e6 << "  ";

      finaloutputfp << g2_counts[runcount] << "  ";

      /*g2*/
      finaloutputfp << setprecision(8) << g2_counts[runcount] / average_n
		    << "  ";
      /*g2 error*/
      finaloutputfp << setprecision(8)
		    << sqrt((double) g2_counts[runcount]) /  average_n << "\n";
    }

    delete [] g2_times;

    delete [] g2_counts;

    finaloutputfp.close();

    cout << "the g2's and errors have been calculated and written to a file"
	 << endl;

    cout << "completed data" << file_counter << endl << endl;
  }

  cout << "completed all files" << endl;

  return 0;

}

void press_enter_continue()
{
  cout << "Press ENTER to continue... " << flush;

	cin.ignore(numeric_limits <streamsize>::max(), '\n');
}
