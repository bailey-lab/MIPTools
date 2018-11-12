

Installation should be fairly easy--the html help document has general install instructions towards the end.  I definitely recommend UNIX over running it on Windows.  If you have no choice however, I have manage to iron out most of the Windows inadequacies/bugs. 


To install the program and the example script:

1) Install Perl if not installed (most standard unix installs have all that is required).
a. For unix see www.perl.org, although for most unix computers Perl is part of the standard install. 
b. For MS Windows, try ActiveState Perl (http://www.activestate.com/Products/ActivePerl/ ). They have already compiled a binary so you don’t.
2) Uncompress the file parasight_v7.4.zip file creating the directory parasight_v7.4 . 
3) Directly within the parasight_v7.4 directory is the parasight executable, it should either be moved or linked into your bin path for ease of use.  On windows machines that have installed ActiveState the simplest solution is to place it in the Perl bin directory (usually C:\Perl\bin).   Parasight code no longer requires additional modules except for Tk.
4) Check to see if parasight runs by typing parasight at a command line.  You should get a summary of options.
5) If it doesn’t work you may need to fix the path or install any modules such as Tk but only if it complains that they are not found.
6) Once you get parasight to run  (i.e. list it’s main options when run without any arguments), try changing to the examples directory and running parasight_examples1.pl –example 1   This is a cheesy scripted tutorial that will demonstrate some of the things parasight can do. This program won’t run unless it can find parasight (i.e. after you put it in a bin directory).  You need to execute parasight_examples1.pl in its directory so it can find the example data.  
7) The other examples aren’t really scripted but they are examples to give you ideas of what parasight can do and how to go about getting parasight doing it.


Enjoy, Jeff

jab@cwru.edu

