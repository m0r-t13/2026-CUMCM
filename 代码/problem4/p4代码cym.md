2026 CUMCM B题 问题4 机器狗自动搜索/定位/清除程序说明
对应文件：problem4_robot_optimized_v2_200s.py

1. 程序概述
本程序用于 2026 年 CUMCM B 题问题 4：目标区域同时存在全向干扰源与定向干扰源，其中定向源有效覆盖角为 180°，方向未知。程序控制机器狗完成：

自动搜索未知频道；

对发现信号进行测向；

利用方向信息定位干扰源；

在满足清除条件后执行 /clear；

对未发现频道在完成保证搜索后判空；

记录 JSONL 日志与 summary 统计结果。

程序全程采用串行 HTTP 请求，兼容全向源，并针对定向源设计了专门的搜索与定位策略。

2. 核心策略
2.1 25 点三角剖分保证搜索网
搜索阶段采用“中心 + 12 个内环点 + 12 个外环点”的 25 点搜索网：

中心点：(0, 0)

内环：12 个点，半径 950 m

外环：12 个点，半径为

R
outer
=
1800
cos
⁡
15
∘
R 
outer
​
 = 
cos15 
∘
 
1800
​
 
即外切目标圆的正 12 边形顶点。

该搜索网的关键性质：

所有剖分三角形边长均小于 1000 m；

任一干扰源都位于某个三角形内，且三个顶点到源的距离均小于 1000 m；

对 180° 定向覆盖半平面，任意经过源的半平面至少包含一个搜索点；

因此，完整扫描 25 点后，定向源必定至少被发现一次。

搜索路线先生成中心点，其余 24 点采用最近邻 + 2-opt 开放路径优化，减少移动距离。

2.2 搜索阶段频道调度
搜索阶段只扫描尚未发现的频道。完成 25 点后，仍未发现的频道才判定为空频道。

程序采用“频道夹心调度”：

若当前频道仍未知，先在新站测它，避免切台；

中间对已知源做零移动同点顺测；

最后扫描剩余未知频道，使站点结束时尽量仍停留在未知频道。

该策略减少频道切换次数，同时不减少保证搜索测量。

2.3 定向源定位与清除
首次发现源后保存示向信息。定位阶段不再假设“1000 m 内必有信号”，而采用方向感知策略：

从最近一次成功测向点沿示向方向前进；

同时向左右偏置形成对称探测点；

若两侧都无信号，则按 1/2 缩小步长；

缩步直到重新进入定向覆盖半平面。

每次成功测向仍使用问题 1 的 ±1° 角楔交会更新保守可行域：

可行域初始为目标区域外接正 180 边形；

用 clip_wedge 保留测得方位附近的角楔；

用 min_forward 约束真实距离大于近距区；

用 MAX_RECEIVE_R = 1500 m 外接多边形保守约束。

当最小包围圆 MEC 半径 <= 19 m 时，直接尝试 /clear；若探测进入 5 m 近距区，则直接清除。

2.4 局部保证清除网格
对于小可行域，程序不一定继续测向，而是构造“定向包围盒网格”：

取可行域凸多边形直径端点作为主轴；

构造旋转包围矩形；

网格单元边长不超过 27.5 m；

单元中心到单元内任一点最大距离约为

27.5
/
2
≈
19.445
m
<
20
m
27.5/ 
2
​
 ≈19.445m<20m
因此每个网格中心点作为 /clear 点可保证覆盖。

若所需网格点数不超过 8，程序直接执行局部清除；否则继续定向测量以缩小可行域。

2.5 同点顺测与全局任务规划
为减少移动成本，程序加入：

搜索站同点顺测：不增加移动，只在预计收益足够大时补测已发现源，每站最多 4 个频道；

顺路清除：若某源已 READY，且绕行额外距离不超过 415 m，则在搜索途中清除；

全局 2-opt 任务规划：在定位阶段，对所有待处理源的下一实际动作点做开放 2-opt，减少多源之间折返。

3. 主要参数
参数	默认值	说明
ARENA_R	1800.0	目标区域半径
MIN_RECEIVE_R	1000.0	最小可接收距离
MAX_RECEIVE_R	1500.0	最大接收距离
BEARING_ERR_DEG	1.0	示向角误差
NEAR_R	5.0	近距直接清除距离
CLEAR_R	20.0	清除半径
READY_RADIUS	19.0	MEC 可清除半径，留 1 m 余量
CIRCLE_SIDES	180	用外接正多边形近似圆盘
P4_INNER_R	950.0	问题 4 内环半径
P4_OUTER_R	1800 / cos(15°)	问题 4 外环半径
P4_RING_N	12	内外环点数
P4_PROBE_FORWARD	240.0	定向追踪最大前进距离
P4_PROBE_LATERAL	30.0	定向追踪最大侧偏距离
P4_MIN_FORWARD	40.0	最小前进距离
P4_MIN_LATERAL	12.0	最小侧偏距离
P4_OPPORTUNISTIC_MAX_PER_STATION	4	每搜索站最多同点顺测频道数
P4_CLEAR_GRID_STEP	27.5	局部清除网格步长
P4_CLEAR_GRID_MAX_POINTS	8	网格清除最大点数
P4_INLINE_CLEAR_MAX_DETOUR_M	415.0	顺路清除最大额外绕行
MOVE_SPEED_MPS	5.0	机器狗移动速度
P4_MIN_PROBE_SCALE	1/64	定向探测最小缩步比例
P4_MAX_SUCCESSFUL_PROBES	12	单源最大成功探测次数告警阈值
P4_MAX_TOTAL_PROBES_PER_SOURCE	60	单源最大总探测次数
MAX_REAL_MARGIN_S	20.0	现实时间安全余量
4. 代码结构
模块/类/函数	作用
基础几何	向量加减、点积、叉积、距离、质心、方位角、半平面裁剪、角楔裁剪、外接多边形等
minimum_enclosing_circle	随机增量法求最小包围圆 MEC
Observation	单次观测记录：位置与示向角
SourceState	单个频道的可行域、观测、清除状态
p4_guarantee_search_points	生成 25 点保证搜索网
p4_search_route	生成 25 点搜索路线
p4_mesh_max_edge	计算三角网最大剖分边长
p4_discrete_directional_gap	离散验证 180° 定向覆盖保证
SimulatorClient	HTTP 客户端，封装 /enter、/measure、/clear、/exit
Problem4Controller	问题 4 主控制器，负责搜索、定位、清除、统计
self_test	本地几何自检
main	命令行入口、日志与 summary 输出
5. 运行方法
5.1 本地自检
不连接模拟器，仅检查几何逻辑：

bash
python problem4_robot_optimized_v2_200s.py --self-test
5.2 正式运行
bash
python problem4_robot_optimized_v2_200s.py --robot-id <参赛队号>
可指定模拟器地址与日志目录：

bash
python problem4_robot_optimized_v2_200s.py \
  --robot-id <参赛队号> \
  --base-url http://127.0.0.1:2026 \
  --log-dir logs
也可通过环境变量设置：

bash
export ROBOT_ID=<参赛队号>
export SIM_BASE_URL=http://127.0.0.1:2026
5.3 运行前注意
请先在模拟器中登录，并启动【问题 4 演练测试/正式测试】，等待界面显示接口已就绪后，再运行本程序。

6. HTTP 接口与重试
程序使用以下接口：

接口	作用
/enter	进入模拟器
/measure	在指定位置测量指定频道
/clear	在指定位置清除指定频道
/exit	退出模拟器
网络请求特点：

全程串行 HTTP；

失败重试时复用同一 request_id；

对 400/404/405/409/413/415 等参数类错误不盲目重试；

对网络错误、超时等进行指数退避重试；

所有请求与响应写入 JSONL 日志。

7. 输出与日志
运行后会在日志目录生成：

text
logs/problem4_YYYYMMDD_HHMMSS.jsonl
logs/problem4_YYYYMMDD_HHMMSS_summary.json
summary.json 主要字段包括：

problem：问题编号；

found_channels / found_count：发现频道；

cleared_channels / cleared_count：清除频道；

clear_rate_vs_found：清除率；

empty_channels：判空频道；

search_strategy：搜索策略统计；

localize_strategy：定位策略统计；

route_metrics：移动距离、估计移动时间、测量次数、切台次数、清除次数等；

virtual_time_s：模拟器虚拟时间；

local_real_elapsed_s：本地现实耗时；

source_details：每个源的观测次数、直径、MEC 半径、清除尝试次数等。

程序返回码：

若 cleared_count == found_count 且 found_count >= 10，返回 0；

否则返回 1，便于演练识别未完全清除情况。

8. 本地自检内容
--self-test 会检查：

25 点保证搜索点数量是否为 25；

最大剖分边长是否严格小于 1000 m；

离散方向最大角缺口是否不超过 180°；

固定搜索路线长度是否可计算；

角域交会后的可行域是否非空、直径有限；

定向追踪满尺度单次试探位移是否小于 1000 m。

自检通过后输出：

text
自检通过。
9. 设计特点总结
保证性搜索：25 点三角网严格覆盖目标圆，定向源 180° 覆盖角不会漏检；

保守定位：使用角楔、外接多边形、最小前向距离等凸约束，避免可行域错误收缩；

方向感知探测：左右对称探测 + 失败缩步，适应定向源半平面覆盖特性；

低移动优化：同点顺测、顺路清除、全局 2-opt 任务规划；

局部保证清除：小可行域用旋转网格覆盖，减少无效测向；

鲁棒通信：重试复用 request_id，JSONL 全量日志；

时间保护：现实剩余时间不足时主动停止新动作，避免超时。

