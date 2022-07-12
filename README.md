

# LatticeGraph.jl

LatticeGraph.jl是一个通过图的方法来构造紧束缚模型的julia包。

## 安装

```julia
pkg> add https://github.com/rainshed/LatticeGraph.jl
```



## 原理

对于一个紧束缚模型，我们可以用一个图（LatGraph）来表示它。在这个图中的每一个点（Node）表示一个自由度（比如单带模型中的一个site），每一条边（Edge）表示两个自由度之间的hopping强度。这个图的连接矩阵就是该模型的哈密顿量的矩阵。通过这种方法，我们可以直观地构造出任意维度任意边界条件的紧束缚模型。

例如，最简单的左右hopping不一样的非厄米哈密顿量$H = \sum (rc^\dagger_ic_{i+1}+sc^\dagger_{i+1}c_i)$的图如下。

![image-20220712171341020](/Users/yupengwang/Library/CloudStorage/OneDrive-shanghaitech.edu.cn/learn/project/NH/LatticeGraph/README.assets/image-20220712171341020.png)

## 生成一张图

### 一维链

下面以长度为100的一维链$H = \sum (rc^\dagger_ic_{i+1}+sc^\dagger_{i+1}c_i)$为例，给出生成一张图的流程。

首先调用LatticeGraph

```julia
using LatticeGraph
```

然后使用函数`LatGraph(f,size)`初始化一个只有Node没有Edge的图，其中`size`代表这个图的最大尺寸，如size = [100,100]代表这个图在100*100的正方形中，每一个整数坐标如(25,79)都可以放置一个Node。`f`是个挑选出你想要的坐标放置Node的函数。比如在我们这个问题中，我们想要一个长度为5的一维链，就可以

```julia
julia> g1 = LatGraph(x->true,[5])
one dimension chain with 5 sites
```

当我们想要在最大尺寸`size`中的每一个点都放上一个Node时，也可以省略第一个参数，即

```julia
julia> g2 = LatGraph([5])
one dimension chain with 5 sites
```

然后调用`AddEdge!(graph,relapos,weight,BC)`函数来给我们的图添上Edge，这个函数会给`graph`中所有相对位置为`repapos`的两个Node之间加上一条权重为`weight`的**单向边**，边界条件`BC`可以在不同方向（当维度大于2的时候）选择周期性边界条件`"PBC"`或者开边界条件`"OBC"`。在这里我们让`g1`为周期性边界条件，`g2`为开边界条件。

```julia
julia> s=1.2;r=0.8;
julia> AddEdge!(g1,[1],r,["PBC"]);
julia> AddEdge!(g1,[-1],s,["PBC"]);
julia> AddEdge!(g2,[-1],s,["OBC"]);
julia> AddEdge!(g2,[1],r,["OBC"]);
```

此时，我们就已经生成了一维链$H = \sum (rc^\dagger_ic_{i+1}+sc^\dagger_{i+1}c_i)$周期边界条件的LatGraph`g1`和开边界条件的LatGraph`g2`。我们可以通过坐标来查看LatGraph中的Node

```julia
julia> g1[1]
position:[1]
Edges:
->[2]  weight:0.8
->[5]  weight:1.2

julia> g2[1]
position:[1]
Edges:
->[2]  weight:0.8 
```

正如我们所预想的，`g1`中的第一个Node有两条分别指向坐标[2]和坐标[100]的Edge，`g2`则只有一条指向[2]的Edge。

我们直接可以对角化生成出来的LatGraph

```julia
julia> val1,vec1 = eigen(g1);
julia> val2,vec2 = eigen(g2);
```

也可以获得其哈密顿量矩阵

```julia
julia> ConMat(g1)
5×5 Matrix{ComplexF64}:
 0.0+0.0im  0.8+0.0im  0.0+0.0im  0.0+0.0im  1.2+0.0im
 1.2+0.0im  0.0+0.0im  0.8+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.2+0.0im  0.0+0.0im  0.8+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.2+0.0im  0.0+0.0im  0.8+0.0im
 0.8+0.0im  0.0+0.0im  0.0+0.0im  1.2+0.0im  0.0+0.0im

julia> ConMat(g2)
5×5 Matrix{ComplexF64}:
 0.0+0.0im  0.8+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 1.2+0.0im  0.0+0.0im  0.8+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.2+0.0im  0.0+0.0im  0.8+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.2+0.0im  0.0+0.0im  0.8+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  1.2+0.0im  0.0+0.0im
```



### 二维多边形边界条件

`LatGraph(poly,size)`函数还可以用于生成一个任意二维多边形的图，其中`poly`是在大小为`size`的方格中的多边形的顶点组成的向量（沿着边的顺序）。

如

```
julia> LatGraph([[1,1],[1,5],[5,1]],[5,5])
.....
....
...  
..   
.    
```

这就生成了一个顶点为[1,1]，[1,5]和[5,1]的三角形的图。

利用这个函数，我们可以用来研究几何依赖的skin effect。

```julia
function Model(L::Int,shape::String)
    """
        H = 2*cos(kx)+i*sin(ky)
    """
    if shape=="Square"
        G = LatGraph([L,L])
    elseif shape=="LX"
        G = LatGraph([[L÷2+1,1],[1,L÷2+1],[L÷2+1,L],[L,L÷2+1]],[L,L])
    end

    AddEdge!(G,[1,0],1);
    AddEdge!(G,[-1,0],1);
    AddEdge!(G,[0,1],0.5);
    AddEdge!(G,[0,-1],-0.5);

    G
end
```



### 二维多带模型

```julia
function Model(
    t1::Real,
    t2::Real,
    γ::Real,
    N::Integer
    )
    """
    H(k) = 
    
    0        t1+γ*sin(ky)+t2*exp(ikx)
    t1-γ*sin(ky)+t2*exp(-ikx)       0

    OBC on x direction, PBC on y direction
    """ 
    G = LatGraph(x->true,[N,N,2])
    BC = ["OBC","PBC","OBC"]
    # sublattice hopping
    AddEdge!(G,[0,0,1],t1,BC)
    AddEdge!(G,[0,0,-1],t1,BC)
    # x direction hopping
    AddEdge!(G,[1,0,-1],t2,BC)
    AddEdge!(G,[-1,0,1],t2,BC)
    # y direction hopping
    AddEdge!(G,[0,1,1],im*γ,BC)
    AddEdge!(G,[0,-1,-1],im*γ,BC)
    AddEdge!(G,[0,1,-1],-im*γ,BC)
    AddEdge!(G,[0,-1,1],-im*γ,BC)
    G
end
```







