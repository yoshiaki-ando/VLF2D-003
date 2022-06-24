# VLF2D-003

## 説明
* 電離圏は等価導伝率でモデル化(単一周波数)
* 電子密度はIRI2016
* 磁場はIGRF13
* 衝突周波数は簡易的な式
* 大地は表面インピーダンス(単一周波数)
* 送受信点の座標を入力する

## 必要なライブラリ
* [AndoLab](https://github.com/yoshiaki-ando/AndoLab_Library)
* [IRI2016](https://github.com/yoshiaki-ando/IRI2016_Cpp_Wrapper)
* [IGRF13](https://github.com/yoshiaki-ando/IGRF13_cpp)
* Eigen3

## 研究すべき内容
* 磁場を1セル毎に異なった値としている。時間がかかるので、適当な距離で求めて、その間は線形補間などでできないか。
* 等価導伝率テンソルも1セル毎に異なった値としている。間引いて線形補間などはできないか。特に日出、日没時に送受信点で電離圏の変化は現れるか？
* Wait and Spiesの電子密度パラメタ(β、h')で近似するとすれば、どこの電子密度にすれば良いか。
