#include <fstream>
#include "tapsim_reader.h"
#define R_BYTES 104

Result::Result (Results results, unsigned int n) {
  char* fp = results.get_filepath();
  if (results.get_binary()) {
    printf("Binary file.\n");
    ifstream file(fp, ios::binary);
    unsigned int fwd = results.get_head_len() + n * R_BYTES;
    file.seekg(fwd);

    file.read(reinterpret_cast<char*>(&index), sizeof(int));
    file.read(reinterpret_cast<char*>(&id), sizeof(int));
    file.read(reinterpret_cast<char*>(&number), sizeof(int));
    file.read(reinterpret_cast<char*>(&voltage), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<0>(start)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<1>(start)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<2>(start)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<0>(stop)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<1>(stop)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<2>(stop)), sizeof(float));
    file.read(reinterpret_cast<char*>(&tof), sizeof(float));
    file.read(reinterpret_cast<char*>(&probability), sizeof(float));
    file.read(reinterpret_cast<char*>(&potentialBefore), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<0>(fieldBefore)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<1>(fieldBefore)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<2>(fieldBefore)), sizeof(float));
    file.read(reinterpret_cast<char*>(&potentialAfter), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<0>(fieldAfter)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<1>(fieldAfter)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<2>(fieldAfter)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<0>(normal)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<1>(normal)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<2>(normal)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<0>(apex)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<1>(apex)), sizeof(float));
    file.read(reinterpret_cast<char*>(&get<2>(apex)), sizeof(float));

    file.close();
  }
  else {
    ifstream file(fp);
    string line;
    unsigned int beg = 0;
    unsigned int end = 0;
    unsigned int fwd = results.get_head_lines();
    for (int i = 0; i <= fwd + n; i++) {
      getline(file, line);
    }

    getline(file, line);

    end = line.find("\t", beg);
    index = stoi(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    id = stoi(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    number = stoul(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    voltage = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<0>(start) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<1>(start) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<2>(start) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<0>(stop) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<1>(stop) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<2>(stop) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    tof = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    probability = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    potentialBefore = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<0>(fieldBefore) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<1>(fieldBefore) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<2>(fieldBefore) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    potentialAfter = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<0>(fieldAfter) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<1>(fieldAfter) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<2>(fieldAfter) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<0>(normal) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<1>(normal) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<2>(normal) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<0>(apex) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<1>(apex) = stof(line.substr(beg, end));
    beg = end + 1;
    end = line.find("\t", beg);
    get<2>(apex) = stof(line.substr(beg, end));
    beg = end + 1;

  file.close();
  }
}
void Result::print(void) {
  printf("Index\tID\tNumber\n");
  printf("%i\t%i\t%i\n", index, id, number);
  printf("Voltage\n");
  printf("%e\n", voltage);
  printf("Start X\t\tStart Y\t\tStart Z\n");
  printf("%e\t%e\t%e\n", get<0>(start), get<1>(start), get<2>(start));
  printf("Stop X\t\tStop Y\t\tStop Z\n");
  printf("%e\t%e\t%e\n", get<0>(stop), get<1>(stop), get<2>(stop));
  printf("TOF\t\tProbability\tPotential Before\n");
  printf("%e\t%e\t%e\n", tof, probability, potentialBefore);
  printf("Field Before X\tField Before Y\tField Before Z\n");
  printf("%e\t%e\t%e\n", get<0>(fieldBefore), get<1>(fieldBefore),
         get<2>(fieldBefore));
  printf("Potential After\n");
  printf("%e\n", potentialAfter);
  printf("Field After X\tField After Y\tField After Z\n");
  printf("%e\t%e\t%e\n", get<0>(fieldAfter), get<1>(fieldAfter),
         get<2>(fieldAfter));
  printf("Normal X\tNormal Y\tNormal Z\n");
  printf("%e\t%e\t%e\n", get<0>(normal), get<1>(normal), get<2>(normal));
  printf("Apex X\t\tApex Y\t\tApex Z\n");
  printf("%e\t%e\t%e\n", get<0>(apex), get<1>(apex), get<2>(apex));
}

Tapsim::Tapsim (char* fp) {
  filepath = fp;
  ifstream file(fp);
  head_len = 0;
  head_lines = 0;
  string buff;
  unsigned int pos = 0;
  bool head = true;
  int i = 0;
  for (i = 0; head && i < 50; i++) {
    getline(file, buff);
    if (buff.find("ASCII") != string::npos) {
      head = false;
      binary = false;
      break;
    }
    else if (buff.find("BINARY") != string::npos) {
      head = false;
      binary = true;
      break;
    }
    else {
      head_len += buff.length();
    }
  }
  head_lines = i;
  if (i >= 50) {
    printf("Error reading header.\n");
    printf("Head Length Read: %u\n", head_len);
    return;
  }
  if (binary == true) {
    string s_entries = buff.substr(7);
    entries = stoi(s_entries);
  }
  else {
    for (i = 0; !file.eof(); i++) {
      getline(file, buff);
    }
    entries = i - 1;
  }
}
void Tapsim::print() {
  printf("Filepath: %s\n", filepath);
  printf("Header Length: %u\n", head_len);
  printf("Header Lines: %u\n", head_lines);
  printf("Binary: %s\n", binary ? "true" : "false");
  printf("Entries: %u\n", entries);
}

Results::Results (char* fp) : Tapsim (fp) {
  for (int i = 0; i < this->get_entries(); i++) {
    results.push_back(Result(*this, i));
  }
}
void Results::print_result(unsigned int n) {
  Result result = this->get_result(n);
  result.print();
}

DataFrame readResults(char* fp) {
  Results c_data(fp);
  StringVector col_names(26);

  DataFrame r_data = DataFrame::create(Named("Index") = index_vec,
                                       Named("ID") =  id_vec,
                                       Named("Number") = number_vec,
                                       Named("Voltage") = voltage_vec,
                                       Named("StartX") = startX_vec,
                                       Named("StartY") = startY_vec,
                                       Named("StartZ") = startZ_vec,
                                       Named("StopX") = stopX_vec,
                                       Named("StopY") = stopY_vec,
                                       Named("StopZ") = stopZ_vec,
                                       Named("TOF") = tof_vec,
                                       Named("Prob") = prob_vec,
                                       Named("PotBef") = potentialBefore_vec,
                                       Named("FieldBefX") = fieldBeforeX_vec,
                                       Named("FieldBefY") = fieldBeforeY_vec,
                                       Named("FieldBefZ") = fieldBeforeZ_vec,
                                       Named("PotAft") = potentialAfter_vec,
                                       Named("FieldAftX") = fieldAfterX_vec,
                                       Named("FieldAftY") = fieldAfterY_vec,
                                       Named("FieldAftZ") = fieldAfterZ_vec,
                                       Named("NormalX") = normalX_vec,
                                       Named("NormalY") = normalY_vec,
                                       Named("NormalZ") = normalZ_vec,
                                       Named("ApexX") = apexX_vec,
                                       Named("ApexY") = apexY_vec,
                                       Named("ApexZ") = apexZ_vec
                                       );
}

// [[Rcpp::export]]
DataFrame readTAPSim(char* fp, string type) {
  if (type == "results") {
    DataFrame results = readResults(fp);
    return results;
  }
  else {
    Rprintf("Not yet implemented.");
  }
}
