import sys

myword = sys.argv[1]

print(myword)

with open('test_%s.txt' % myword, 'w') as f:
    f.write('Testing: %s' % myword)