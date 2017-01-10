void ConverHexToInt(){
	//std::string s = "0xfffefffe";
        std::string s = "cccc00";
    unsigned int x = std::stoul(s, nullptr, 16);
    cout << "\n"<<x<<"\n";
}
