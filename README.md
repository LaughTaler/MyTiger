<<<<<<< HEAD
# TiGER 核心库


## 规范


### 注释规范

注释使用Doxgen统一管理：https://www.doxygen.nl/index.html

从源代码生成文档
Doxygen是事实上的标准工具，用于从注释C++源生成文档，但它也支持其他流行的编程语言，如C，Objective-C，C#，PHP，Java，Python，IDL（Corba，Microsoft和UNO / OpenOffice风格），Fortran，以及在某种程度上的D.Doxygen还支持硬件描述语言VHDL。

Doxygen可以通过三种方式为您提供帮助：

它可以从一组记录的源文件生成在线文档浏览器（HTML格式）和/或离线参考手册（在$\mbox{\LaTeX}$中）。还支持在RTF（MS-Word），PostScript，超链接PDF，压缩HTML和Unix手册页中生成输出。文档直接从源代码中提取，这使得保持文档与源代码一致变得更加容易。
您可以配置 doxygen 以从未记录的源文件中提取代码结构。这对于在大型源代码发行版中快速找到自己的方式非常有用。Doxygen还可以通过包含依赖图，继承图和协作图来可视化各种元素之间的关系，这些图都是自动生成的。
您也可以使用 doxygen 创建普通文档（就像我为 doxygen 用户手册和网站所做的那样）。
Doxygen是在Mac OS X和Linux下开发的，但设置为高度可移植。因此，它也运行在大多数其他Unix风格上。此外，Windows的可执行文件也可用。


### 代码规范

命名规范参考Google规范手册：

https://google-styleguide.readthedocs.io/zh_CN/latest/google-cpp-styleguide/contents.html

C++ 是 Google 大部分开源项目的主要编程语言. 正如每个 C++ 程序员都知道的, C++ 有很多强大的特性, 但这种强大不可避免的导致它走向复杂，使代码更容易产生 bug, 难以阅读和维护.

本指南的目的是通过详细阐述 C++ 注意事项来驾驭其复杂性. 这些规则在保证代码易于管理的同时, 也能高效使用 C++ 的语言特性.

风格, 亦被称作可读性, 也就是指导 C++ 编程的约定. 使用术语 “风格” 有些用词不当, 因为这些习惯远不止源代码文件格式化这么简单.

使代码易于管理的方法之一是加强代码一致性. 让任何程序员都可以快速读懂你的代码这点非常重要. 保持统一编程风格并遵守约定意味着可以很容易根据 “模式匹配” 规则来推断各种标识符的含义. 创建通用, 必需的习惯用语和模式可以使代码更容易理解. 在一些情况下可能有充分的理由改变某些编程风格, 但我们还是应该遵循一致性原则，尽量不这么做.

本指南的另一个观点是 C++ 特性的臃肿. C++ 是一门包含大量高级特性的庞大语言. 某些情况下, 我们会限制甚至禁止使用某些特性. 这么做是为了保持代码清爽, 避免这些特性可能导致的各种问题. 指南中列举了这类特性, 并解释为什么这些特性被限制使用.

### TiGER 接口管理规范

* 接口从上到下顺序为输入参数/输出参数
* 输入参数只能使用const引用类型
* 输出参数只能使用引用类型
* 不允许使用任何形式的指针,指针应当使用引用替代
* 尺寸函数使用std::function
* 所有选项参数使用std::string传导，其选项使用CLI解析
* 所有非固定长度数组类使用std::vector
* 所有固定长度数组类使用std::array
* 不允许使用非上述类型的标准库类型
* 不允许使用基本数据类型，基本数据类型应当视为选项参数
* 必须在C++11国产机器上能编译

=======
# MyTiger
a repository for update and get code of tiger
>>>>>>> 3e5f6c4eb43b8b1f41c917b6503cf36d68a6b66d
