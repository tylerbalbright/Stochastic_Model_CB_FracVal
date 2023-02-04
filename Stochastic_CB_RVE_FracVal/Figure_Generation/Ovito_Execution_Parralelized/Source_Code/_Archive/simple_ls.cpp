#include <string>
#include <iostream>
#include <boost/filesystem.hpp>
#include <unistd.h>
#include <vector>

int main()
{
    std::vector<std::string> file_names;
    std::string cwd = getcwd(NULL,0);
    std::string dir = cwd + "/_Archive/";
    boost::filesystem::path p(dir);
    for (auto i = boost::filesystem::directory_iterator(p); i != boost::filesystem::directory_iterator(); i++)
    {
        if (!boost::filesystem::is_directory(i->path())) //we eliminate directories in a list
        {
	    file_names.push_back(i->path().filename().string());
            std::cout << i->path().filename().string() << std::endl;
        }
        else
	{
            continue;
    	}	
    }	
}
