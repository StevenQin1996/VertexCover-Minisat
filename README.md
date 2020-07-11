# ECE650 : Final Project

Replace FIRST_NAME, LAST_NAME, WATIAM, and STUDENT_ID, EMAIL in
`user.yml` by your first and last name, WatIAM id, student number, and
email, respectively. The information must be entered for all members
of your tema. Whenever possible, please use ASCII characters.

Do not change the format `user.yml`. We will parse it
automatically. Only enter the information requested.

The main executable file for your solution to the assignment should be
`ece650-prj`.

Do not include MiniSat in your repository. We will clone it in your
repository using the command. This is exactly the same as was done in
Assignment 4.

```
git clone https://github.com/agurfinkel/minisat
```

Note that for the project you have to create a `CMakeLists.txt` on
your own. You can use examples from previous assignments or from
course examples on GitHub.

Do not forget to include your report in `report.pdf`

Commit your changes and submit on GitHub.

我尝试了在头文件里加上三个private的vector 分别用来记录三个算法的结果。在头文件里的public function 来set和get这三个VECTOR的值。但是发现工作量不小。需要把我们的现在已知所有的function全部转移到头文件里才能简洁的运行。不这么做的话只能把代码强拼在一起，之后debug会很麻烦。

我觉得还是试试看通过在函数中通过指针来修改reference的值，但似乎需要在本地建立三个vector来记录三个算法的结果。每个算法返回一个vector。最后在thread1中输出这三个ector的值。

还有 就是我觉得我们需要用互斥锁来保证四个线程的运行顺序。理想的运行顺序应该是1(input)->2->3->4->1(output).
