#ifndef SIMSUITE_TESTS_CORETESTS_HPP_
#define SIMSUITE_TESTS_CORETESTS_HPP_

#include <string>

void TEST_isomorphicNewick_TRUE(void);
void TEST_ALL(void);
bool ASSERT(bool assertion, std::string title, std::string info);

#endif