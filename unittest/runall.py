import unittest
loader = unittest.TestLoader()
suites = loader.discover('.', pattern='test_*.py')
result = unittest.TestResult()
for suite in suites:
    suite.run(result)

for error in result.errors:
    print error[0]
    print error[1]
    print '\n'

for failure in result.failures:
    print failure[0]
    print failure[1]
    print '\n'

print result


