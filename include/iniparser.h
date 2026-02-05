/*
this header contain function prototypes that helps to parse a .ini file using STL map container and template functions 
*/
#ifndef INIPARSER_H
#define INIPARSER_H

#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <species.h>
#include <tuple>

class Domain;

/**
 * @class INIParser
 * @brief A simple INI file parser that reads configuration files( .ini extension files) and stores the data in a nested map structure.
 */

/// @brief  
class INIParser 
{
    public:
    /**
     * @brief Constructor for the INIParser class.
     * @param domain Reference to the simulation domain.
     */
    INIParser(Domain &domain):domain(domain){};
    /**
     * @brief Parses the given INI file and returns a nested map containing sections and key-value pairs.
     * @param filename Path to the INI file to be parsed.
     * @return A nested map where the outer map's keys are section names and the inner maps contain key-value pairs for that section.
     */
    static std::map<std::string, std::map<std::string, std::string>> parse(const std::string &filename);
    /**
     * @brief Retrieves an integer value from the specified section and key in the parsed INI data.
     * @param section The section name in the INI file.
     * @param key The key name within the section.
     * @return The integer value associated with the specified key in the section.
     */
    static int getInt(const std::map<std::string, std::string> &section, const std::string &key);
    /**
     * @brief Retrieves a double value from the specified section and key in the parsed INI data.
     * @param section The section name in the INI file.
     * @param key The key name within the section.
     * @return The double value associated with the specified key in the section.
     */
    static double getDouble(const std::map<std::string, std::string> &section, const std::string &key);
    /**
     * @brief Retrieves a string value from the specified section and key in the parsed INI data.
     * @param section The section name in the INI file.
     * @param key The key name within the section.
     * @return The string value associated with the specified key in the section.
     */
    static std::string getString(const std::map<std::string, std::string> &section, const std::string &key);
    /**
     * @brief Splits a string into a vector of substrings based on the specified delimiter.
     * @param str The string to be split.
     * @param delimiter The character used to split the string.
     * @return A vector of substrings.
     */
    static std::vector<std::string> split(const std::string &str, char delimiter);
    /**
     * @brief Extracts load type to initialize particle initial pistion and velocity.
     * @param position_init The position initialization string (e.g., "uniform , random, 0.5sin(10) etc).
     * @return A tuple containing the load type and two  extra parameter.
     */
    static std::tuple<std::string, int, double> loadtypeextract(const std::string &position_init);
    /**
     * @brief Parses a collision group string into a vector of integer pairs representing collision groups.
     * @param collGroupStr The collision group string (e.g., "0-1,2-3").
     * @return A vector of pairs of integers, where each pair represents a collision group.
     */
    static std::vector<std::pair<int, int>> parseCollGroup(const std::string &collGroupStr);
    
    private:
    /// @brief Helper function to trim whitespace from both ends of a string.
    static void trim(std::string& str);
    /// @brief Reference to the simulation domain.
    Domain &domain;
    
};

#endif // INIPARSER_H
