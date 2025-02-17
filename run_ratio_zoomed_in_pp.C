#include "createratio.C"

void extractRangesAndCallFunction(const std::string &filename, const std::string &rootFilename)
{
  std::ifstream file(filename);
  if (!file)
  {
    std::cerr << "Error: Unable to open file " << filename << std::endl;
    return;
  }

  std::vector<double> ranges(4); // Stores xmin, xmax, ymin, ymax
  std::string line;

  while (std::getline(file, line))
  {
    std::istringstream iss(line);
    std::string key;
    if (line.find("X-axis Range:") != std::string::npos)
    {
      char ch;
      iss >> key >> key >> ch >> ranges[0] >> ch >> ranges[1]; // Extracts xmin, xmax
    }
    else if (line.find("Y-axis Range:") != std::string::npos)
    {
      char ch;
      iss >> key >> key >> ch >> ranges[2] >> ch >> ranges[3]; // Extracts ymin, ymax
    }
  }

  file.close();

  // Debug print
  std::cout << "Extracted ranges: [" << ranges[0] << ", " << ranges[1] << ", " << ranges[2] << ", " << ranges[3] << "]" << std::endl;

  // Call createratio with the extracted values
  createratio(rootFilename, ranges[0], ranges[1], ranges[2], ranges[3]);
}

void run_ratio_zoomed_in_pp()
{
  // format: Name, shiftlow, shifthigh, smearlow, smearhigh

  // I only need to zoom in PbPb, pp inclusive.

  TString variation[7] = {"nominal", "tnpU", "tnpD", "acooff", "nominal_no_bk_sub", "nominal_binning", "nominal_range"};

  for (int var = 0; var < 7; var++)
  {
    for (int iseta = 0; iseta < 2; iseta++)
    {
      TString filename = "";
      TString type = "";
      TString savedname = "";

      if (iseta == 1)
        type = "eta_";
      if (iseta == 0)
        type = "FA_";

      filename = TString("../pp_18/zoomin/") + TString("pp_") + type + variation[var] + "_" + "22.txt";
      savedname = TString("./zoomin_pp/") + TString("pp_") + type + variation[var] + "_" + "22.root";

      extractRangesAndCallFunction(std::string(filename.Data()), std::string(savedname.Data()));
    }
  }

  cout << "Finished" << endl;
}