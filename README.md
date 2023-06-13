## 简介
- 分享一些我在完成本科毕业设计的过程中用到的代码
- 包括非常基础的R语言代码和linux (zsh shell script)
- 用这些代码实现了以下功能
    - 差异表达分析（edgeR）
    - 热图、火山图绘制（pheatmap、EnhancedVolcano）
    - GO与KEGG富集分析及气泡图绘制（clusterprofiler）
    - 文件整理（linux、R共同实现）
- 差异表达分析和富集分析属于是生信上古遗物了，没什么新意，毕竟本科生搞的东西，主要是帮自己熟悉一下R语言和linux
- 唯一比较特殊的一点是我用的生物没有现成的注释信息，所以KEGG和GO库都是自己整理出来的自建库。

## 代码文件信息
- `DiffExpAnalysis.R`：差异表达分析、差异基因的火山图、热图绘制。
- `KEGGArrange.sh`：整理KEGG自建库的文件（GO类似，不放了），本质上是处理文件的linux脚本，没什么逻辑，而且不能直接跑（因为我当时是写好代码后从里面一段段复制粘贴到终端跑的，懒得整理成完整的脚本了）。每一步都有简单的注释，主要目的是帮自己熟悉一下linux代码。
- `kegg111.R`：整理文件的过程中有一步用了R，在`KEGGArrange.sh`的注释中有说明。很烂的脚本，用R语言写for循环真的很搞笑。
- `KEGG0524.R`：富集分析和气泡图绘制，用到的文件`ko_fin.txt`和`ko_anno.txt`的格式示例放在这里了。这个代码KEGG和GO或者其他数据库的富集分析都是通用的，物种是非模式生物，用的是自建库。
