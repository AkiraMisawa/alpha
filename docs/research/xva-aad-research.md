# AADによるXVAリスク計算アーキテクチャ 調査レポート

調査日: 2026-11-13。Matt Pocock の research スキル（mattpocock/skills, `skills/engineering/research`）の手順に従い、一次資料を調べて出典付きでまとめた。第1部は統合した結論、第2部と第3部は調査ノートの全文（出典はそこに記載）。

表記: [確認] は原典で確かめた内容、[抄録のみ] は本文未読、[推測] は本レポートの推論、[未確認] は原典に当たれなかったもの。【ベンダー公表値】は第三者の検証がない数値。

リポジトリ AkiraMisawa/alpha は調査時点でコミットのない空の状態だったため、既存コードとの照合はしていない。

---

# 第1部 結論と設計への示唆

## 1. 合意済み方針の判定

| 合意事項 | 判定 | 根拠と注意 |
|---|---|---|
| 回帰と集計をネッティングセット単位の1タスクに統合 | 支持 | 回帰の随伴は全パスにまたがるリダクション（Ψ̄, Ω̄）を要するので、全パスが同じタスクにある構成と整合する（Antonov 2016、Capriotti–Jiang–Macrina 2017）。ネッティングセットの大きさの偏りがタスク負荷の偏りになる。 |
| 集約ツリーと、非線形集約でのベクトル受け渡し・シード配布 | 支持 | KVAの資本はポートフォリオ水準の要件（Green–Kenyon–Dennis 2014）、SIMMは平方根集約（ISDA SIMM）。スカラー合算では上位の値にならない。ツリーの根にメモリと負荷が集中する点に注意。 |
| **逆伝播でβを定数として固定** | **要修正** | Capriotti–Jiang–Macrina (2017) は、回帰値をそのまま露出に使うXVAではβ固定は「まったく使えない」とし、ボラ感応度がバンプと大きく乖離したと報告。固定してよいのは、回帰を最適行使の判定だけに使う場合の1次感応度に限る（Fries 2019）。 |
| チェックポイント | 支持 | Capriotti & Giles 2012、Savine 2018 第13章、Griewank & Walther の revolve、NAG の事例（メモリ不足→8.7MB）。 |
| ハイブリッドモデル、線形商品は解析式、回帰はエキゾチックのみ | 支持 | 解析式とAADは矛盾しない（Capriotti & Lee 2014）。回帰の範囲を絞ることは、β問題の影響範囲も絞る。 |
| キャリブレーションを陰関数定理で微分 | 支持 | Henrard 2014、Christianson 1994、Goloubentsev–Lakshtanov–Piterbarg 2022。exact-fit でない場合とソルバー収束誤差に注意。 |
| PV調整（adj入力ノード化と固定モード） | 文献なし | 連鎖律としては自明だが、直接の文献は見つからなかった。社内で仕様と検証を文書化する必要がある。 |
| クロスガンマを bump-over-AAD で | 支持 | Capriotti & Giles 2024 が、tangent-over-adjoint と同じコストのオーダーで実装が簡単な方法として挙げる。不連続点の2次感応度は平滑化とセットで扱う。 |
| C++コア + C API + P/Invoke | 支持 | XAD のテープAPIが設計に合う。Unix ではネイティブ例外がC境界を越えられないため、全エクスポート関数で捕捉が必須。 |

## 2. 最重要の論点: βの扱い

βを固定すると、XVAでは回帰した継続価値 H がそのまま露出 V⁺ として入るため、「行使境界が最適に近いので1次の効果が消える」という議論が使えない。選択肢は3つある。

1. **βの随伴を実装する（推奨）。** 前進 β = Ψ⁻¹Ω に対し、逆方向は Ω̄ = Ψ⁻ᵀβ̄、Ψ̄ = −Ω̄βᵀ。基底次元の小さな線形代数で済み、追加コストは小さい（Capriotti–Jiang–Macrina 2017 式(75)）。固定モードはオプションとして残し、差を計測できるようにする。
2. **Fries の方式。** 逆方向で「随伴の条件付き期待値」を取り、回帰そのものは微分しない（Fries 2019 §3.2）。
3. **回帰を判定だけに使う。** 回帰は行使判定と露出の符号判定にだけ使い、露出の値は実現キャッシュフローから組む。この場合は準最適性の議論が使える余地がある [推測。Andreasen 型と言われるが原典未確認]。

どれを選んでも、ボラ感応度（ベガ）をバンプと照合する検証を必須にする。ハザードレート感応度はβの扱いの影響を受けない。

なお ORE の計算グラフでは条件付き期待値の勾配が未実装（`QL_FAIL`）で、OSSの既製実装をそのまま頼ることはできない。

## 3. PV調整の設計案 [本レポートの提案]

定義: adj = F − M(0)、V(t) = M(t) + adj × d(t)。F はフロント時価、M はモデル時価、d は d(0)=1 で満期に向けて0へ減衰する関数。

- **連結モード**: adj をテープ上の入力ノードにする。逆伝播で adj̄ = Σ_{path,t} V̄(path,t) × d(t) が得られ、市場クオート θ への感応度は
  dXVA/dθ = (モデル経由の感応度) + adj̄ × (∂F/∂θ − ∂M(0)/∂θ)。
  ∂F/∂θ はフロントのデルタ、∂M(0)/∂θ は同じテープで得られる。フロントのデルタとモデルの感応度を同じクオート空間に揃えるヤコビアンが要る。
- **固定モード**: adj を定数として扱い、adj̄ を伝播しない。フロントのデルタが手に入らない取引や、切り分け検証に使う。
- 決めるべきこと: adj の粒度（取引単位かネッティングセット単位か）、d(t) の形（線形、想定元本プロファイル比例など）、行使後の扱い、LSM の回帰に adj を入れるかどうか。
- 検証: 連結モードはフロントとモデルを両方バンプした差分と、固定モードは adj を据え置いたバンプと、それぞれ一致することを確かめる。

## 4. ライブラリと実行基盤

推奨順位は XAD → AADC → dco/c++。詳細は第3部。

- **XAD**: `clearAll` / `getPosition` / `resetTo` / `computeAdjointsTo` / nested recording / callback があり、パス単位テープの再利用と、LSM・集約ノードの手書き随伴の受け口が揃う。テープはスレッドローカル。AGPL-3.0 のため商用ライセンスの法務確認が要る。
- **AADC**: 記録1回・再生多数のJITでイントラデイ再計算に強い。商用、分岐の書き換えが必要、GPU出力なし。XADで設計を固めたあとに、ボトルネックで評価する二段構えが妥当。
- **どれを選んでも、LSMと集約ノードはベクトル単位のカスタム随伴として自作する**ことになる。
- **GPU/JAX**: 分岐の少ないデータ並列カーネルで、倍精度（JAXは既定がfloat32）、GPUメモリに収まること、C#からプロセス外呼出を許容できること、が採用条件。日次クロスガンマがCPUグリッドで時間内に終わるなら、CPU AADが保守性とC#統合で優位。

## 4a. 有料ライブラリを使わない場合の方針 [2026-11-13 追記。本レポートの設計判断で、実測比較の出典はない]

前提: 有料ライブラリは使わない。C#は入口だけでよく、コアはC++またはRust。この前提では、第4節の推奨順位（XAD → AADC → dco/c++）は当てはまらない。AADC と dco/c++ は商用で、XAD は AGPL-3.0 のため非公開の業務システムに組み込むと商用ライセンスが必要になることが多い。

**言語: C++ を推奨。**
- 数値型をテンプレートにして同じコードを double とAAD型で動かす定石が素直に書ける。参照実装（Savine の CompFinance、QuantLib、ORE）もC++。
- 無料の選択肢がある: CoDiPack（GPLv3）、CppAD（EPL-2.0 または GPL-2.0+）、Adept（Apache-2.0）、Enzyme（Apache-2.0 with LLVM Exceptions）。
- Rust は実績のあるテープ型AADのクレートが乏しく、Enzyme による自動微分は nightly 限定の実験機能。スレッドごとのテープを共有する書き方も借用規則と相性がよくない。

**テープ: 全パスの配列を1ノードとするテープを自作する。**
- 演算1つがパス数ぶんのSIMDループになり、記録のオーバーヘッドがパス数に比例しない。ORE の計算グラフ（ノードの値がパス数長ベクトル）と同じ発想。
- LSMと集約ノードはどのみち手書きの随伴になるので、この形と自然につながる。
- メモリはノード数 × パス数 × 8バイト。時点チェックポイントと不要ノードの解放が前提。
- C APIの境界は粗い粒度（ネッティングセット1タスク単位）にする。

**進める順序: まず自作テープと基本演算だけで作る。Enzyme は必須ではない。**
1. 自作テープに基本演算（加減乗除、exp、log、max、指示関数など）の随伴を持たせる。その演算で書いた解析式やペイオフはすべて微分できる。外部への依存がなく、MSVCでもビルドできる。
2. プロファイルを取り、遅い箇所だけを1つのループにまとめた関数に置き換える。その随伴は手書きするか、Enzyme で生成する。
- Enzyme を最初から入れると、LLVMのバージョン固定とビルド環境の制約を先に抱える。
- Windows Server での実行: clang-cl でMSVC互換のネイティブDLLを作れば動き、C#からP/Invokeで呼べる。ただし Enzyme 自体のWindowsでの安定動作は [未確認]。Enzyme を使うなら計算ノードはLinuxのほうが安全。

**XAD と Enzyme: 外部ライブラリを本体に入れるなら Enzyme、XAD は照合用。**

| 観点 | XAD | Enzyme |
|---|---|---|
| ライセンス | AGPL-3.0。非公開システムへの組み込みは実質有料 | Apache-2.0 with LLVM Exceptions |
| 方式 | 演算子オーバーロードで実行時にテープを記録 | LLVM IR をコンパイル時に微分 |
| 速度 | テープの記録と再生のコストがかかる | 最適化後のコードを微分するので速くなりやすい [推測] |
| ビルド | 主要コンパイラで使える | Clang/LLVM が前提。MSVC では使えない |
| 制約 | 型の置き換えが必要 | 仮想関数や外部ライブラリ呼び出しで失敗することがある |
| 金融での実績 | QuantLib 統合あり | 未確認 |

- Enzyme の役割は、線形商品の解析式、ペイオフ、キャリブレーション残差など、小さく単純な関数の随伴を作ること。それを自作テープのベクトルノードとして包む。
- XAD や CoDiPack はテストでの照合用に限る。照合用でも社内のライセンス規程の確認は要る。
- 最初の確認: LGMのスワップ評価など小さな関数1つで、自社のビルド環境で Enzyme が通るかを試す。

## 5. メモリとテープ

- 3段構成: (a) 順方向で全パスのパス×時点ベクトルを保持、(b) 回帰・集約はベクトル単位の外部関数、(c) パスごとの局所計算は 記録→逆伝播→巻き戻し（Savine の CompFinance と同じ形）。
- 見積り: ベクトル1本 = パス数 × 時点数 × 8バイト。1万パス×200時点で約16MB、随伴込みで倍。これに上位へ渡す本数を掛ける。テープ側は1パスあたりの記録バイト数（XADなら `getMemory()` で実測）× スレッド数。

## 6. 検証方法

1. Taylor test（残差が O(h²) で減り、収束率が約2になること）。
2. 共通乱数での中心差分との照合。h を対数掃引して誤差の谷を選ぶ。
3. β固定とβ込みの差をバンプと比べて定量化（Capriotti–Jiang–Macrina の再現）。
4. 指示関数の平滑化あり・なしの比較。
5. PV調整の2モードそれぞれの照合（第3節）。
6. 現行システム（バンプ）との並行稼働での突合。

## 7. 段階的な移行計画 [本レポートの提案]

| 段階 | 内容 | 完了条件 |
|---|---|---|
| 0 | 現行バンプの基準データと照合ハーネスを整備（共通乱数、固定シード） | 代表ネッティングセットで現行値を再現 |
| 1 | 線形商品（解析式）のCVAをネッティングセット1タスクでAAD化。パス単位テープ、IFTキャリブレーション | デルタ・ベガがバンプと許容誤差内で一致 |
| 2 | エキゾチックのLSMを同タスクに統合。βの随伴と固定モード、チェックポイント | β込みがバンプと一致、固定モードとの差を記録 |
| 3 | 集約ツリー（ネッティングセット→相手先→ポートフォリオ）、FVA/KVAのベクトル受け渡しとシード配布、PV調整の2モード | グローバルバリアなしで現行XVAを再現 |
| 4 | bump-over-AAD のクロスガンマ、C APIとC#ラッパ、現行との並行稼働 | 日次クロスガンマが時間内に完了 |
| 5 | MVA/SIMM（将来時点のクオート感応度）、イントラデイ化（AADC・GPUの評価、必要なら代理モデル） | イントラデイの目標時間を達成 |

## 8. 未確認の事項

- Savine/Huge「LSM Reloaded」が β をどう扱うか（原典非公開）。
- Antonov 他 (2018) が「回帰の微分を回避できる」とする条件（本文未読）。
- Andreasen "CVA on an iPad Mini" の一次資料。
- XAD と ORE の連携、XAD・dco/c++ の勾配検証ユーティリティの有無。
- bump-over-AAD と forward-over-reverse をXVAで定量比較した原典。
- 書誌の訂正: "Real-time counterparty credit risk management in Monte Carlo" は Capriotti, Lee, Peacock（Risk 2011）。Capriotti, Jiang, Macrina (2017) は "AAD and least-square Monte Carlo: fast Bermudan-style options and XVA Greeks"。

---

# 第2部 文献調査ノート


調査日: 2026-11-13。表記の決まり:
- **[確認]**: 原典、または出版社・著者の公式ページ（抄録、書誌、本文PDF）で確かめた内容。
- **[抄録のみ]**: 抄録や書誌は確かめたが、本文はペイウォールのため読めていないもの。
- **[推測]**: 本ノートの筆者による推論で、文献に書かれた事実ではないもの。
- **[未確認]**: 原典に当たれなかったもの。

---

### 1. Giles & Glasserman と Capriotti 系

**1.1 Giles & Glasserman, "Smoking adjoints: fast Monte Carlo Greeks", *Risk*, 2006年1月** [抄録のみ]
- 工学の設計最適化で使われてきた adjoint 法を、モンテカルロの pathwise 感応度計算に持ち込んだ論文。多数の入力（LMM の初期フォワード曲線など）に対する感応度を、1本の逆方向スイープでまとめて得る。
- AD を解説した付録付きの長い版は、Oxford の技術報告 NA-05-15 にある。
- 出典: https://people.maths.ox.ac.uk/gilesm/codes/libor_AD 、https://www.risk.net/node/1500261
- 適用上の注意: pathwise 法を使うには、ペイオフがパスごとに（ほぼ至る所で）微分可能で、Lipschitz 連続であることが前提になる。

**1.2 Capriotti, "Fast Greeks by algorithmic differentiation", *J. Computational Finance* 14(3), 2011, pp.3–35** [抄録のみ]
- pathwise 法でボトルネックになるペイオフの微分を AD で自動化し、adjoint モードなら従来法より桁違いに速いと主張している。
- 出典: https://www.risk.net/node/2160393
- 不連続ペイオフの扱いは、Capriotti & Giles のレビュー（1.6）が「Capriotti (2011) はペイオフの平滑化を論じている」と紹介している。本文での具体的な扱いは [未確認]。

**1.3 Capriotti & Giles, "Fast correlation Greeks by adjoint algorithmic differentiation", *Risk* 2010（arXiv:1004.1855）** [確認]
- Gaussian copula の全ペア相関に対する感応度を、資産数によらず「価格計算の高々4倍」のコストで計算できると述べている。
- pathwise 法が成り立つにはペイオフが十分に正則（例: Lipschitz）である必要があることも明記している。
- 出典: https://arxiv.org/pdf/1004.1855

**1.4 Capriotti & Giles, "Algorithmic differentiation: adjoint Greeks made easy", *Risk* 2012** [抄録のみ]
- AD を体系的に使って adjoint 計算を実装する手順を示した論文。
- レビュー（1.6）はチェックポイントの参照先としてこの論文を挙げている。
- 出典: https://www.risk.net/node/2201839

**1.5 リアルタイムCCR（依頼文の書誌に食い違いあり）**
- "Real-time counterparty credit risk management in Monte Carlo" の著者は **Capriotti, Lee, Peacock**、掲載は *Risk*、2011年5月25日。AAD で CCR 計算のコストを数百分の一にできると述べている [抄録のみ]。出典: https://www.risk.net/node/2073092
- **Capriotti, Jiang, Macrina** の 2017 年論文は別物で、"AAD and least-square Monte Carlo: fast Bermudan-style options and XVA Greeks", *Algorithmic Finance* 6(1–2), 2017, pp.35–49 である [確認]。出典: https://ideas.repec.org/a/ris/iosalg/0057.html 、本文プレプリント: https://discovery.ucl.ac.uk/id/eprint/1523595/1/Macrina_AADMC_final.pdf （内容は第2節で扱う）
- 同じ3著者には "Real-time risk management: an AAD-PDE approach"（プレプリント 2015年11月）もある。PDE で解く場合に AAD と陰関数定理（IFT）を組み合わせ、キャリブレーション経由の感応度をバンプより約3桁速く計算したと報告している [確認]。掲載誌は [未確認]。出典: https://discovery.ucl.ac.uk/id/eprint/1522613/1/CJM-IJFE2015NOV24.pdf
- Capriotti & Lee, "Adjoint credit risk management", *Risk* 2014: CDS や CDS インデックスのように、半解析モデルで評価する商品にも AAD を適用している [抄録のみ]。出典: https://www.risk.net/node/2356909

**1.6 Capriotti & Giles, "15 Years of Adjoint Algorithmic Differentiation in Finance"（2024年1月22日付レビュー）** [確認]
- 不連続ペイオフについて:
  - pathwise 法の「最大の制約」は、デジタルやバリアのような不連続ペイオフを扱えないことだとしている。
  - 実務的な回避策はヘヴィサイド関数を滑らかな関数に置き換えることだが、感応度に有限のバイアスが入る。
  - 他の選択肢として、Dirac デルタの解析的積分、条件付き期待値の利用、vibrato（尤度比法と組み合わせる方法）を挙げている。
- 2次感応度について:
  - AAD でヘッセ行列を得るコストは入力数に比例する。
  - 「AAD勾配の有限差分」でヘッセ行列を推定するやり方も同じコストのオーダーで、実装が簡単だという利点を挙げている（8節参照）。
- 出典: https://people.maths.ox.ac.uk/gilesm/files/AAD_Review.pdf

**1.7 Capriotti, "Likelihood ratio method and algorithmic differentiation: fast second order Greeks", *Algorithmic Finance* 4(1–2), 2015, pp.81–87** [抄録のみ]
- AAD に pathwise 法と尤度比法（LRM）を組み合わせ、2次グリークスを推定する方法。
- レビューによれば、バスケットのガンマ行列全体をバスケットの大きさによらず価格計算の約4倍で計算している。
- 出典: https://ideas.repec.org/a/ris/iosalg/0038.html

---

### 2. Savine 系と LSM における回帰係数 β の扱い（最重要）

**2.1 Savine, *Modern Computational Finance: AAD and Parallel Simulations*, Wiley 2018** [書誌・著者ページのみ]
- 並列 C++、モンテカルロ、AAD を扱う本。第13章でチェックポイントを扱うことは著者ページで確かめた。
- マルチスレッド AAD（スレッドごとのテープ）は書籍の主題だが、章ごとの記述は本調査では [未確認]。
- 著者ページの紹介文には LSM（回帰）への言及がない。
- 出典: https://antoinesavine.com/books-by-antoine-savine/

**2.2 Savine & Andreasen, *Modern Computational Finance: Scripting for Derivatives and xVA*, Wiley 2021** [目次で確認]
- 目次の主な内容:
  - 第15章 American Monte-Carlo: LSM、proxy、追加の回帰変数、多重行使と再帰。
  - 第IV部（第16–20章）: リスク感応度、平滑化、自動平滑化アルゴリズム、fuzzy logic。
- 目次に xVA 専用の章はなく、xVA が現れるのは書名だけである。
- 「βを固定してよい」という主張が本書にあるかどうかは [未確認]。
- 出典: https://biblioteca.izt.uam.mx/vufind/Record/EBOOKCENTRAL_on1285169849/Details

**2.3 Savine (Huge & Savine), "LSM Reloaded – Differentiate xVA on your iPad Mini"（Oxford セミナー、2017年6月8日）** [抄録のみ]
- 抄録の要点:
  - Andreasen による、それまで完全には公開されていなかった xVA の再定式化を土台にしている。
  - 有担保・無担保の xVA 評価アルゴリズムを示す。
  - シミュレーションと LSM ステップを合わせて微分する際の、メモリ使用量と効率の問題に対処する独自アルゴリズムを提示している。
- 出典: https://www.maths.ox.ac.uk/node/24940
- 回帰係数を固定するのか、微分するのかは抄録には書かれていない [未確認]。

**2.4 Capriotti, Jiang, Macrina (2017)**（書誌は 1.5） [本文で確認]
- 本文の主な主張:
  - バミューダン型では、回帰による行使境界が最適に近いことを根拠に、βを固定する近似が使われる（§3.1）。ただし精度は基底関数に依存する（§3.4）。
  - この近似は「回帰を XVA に用いる場合はまったく使えない」（§3.4）。XVA では回帰で得た保有価値 H がそのまま露出 V⁺/V⁻ として XVA ペイオフに入る（式(43)(44)）ので、準最適性の議論が成り立たない（§4.2）。
  - β を固定した AAD で求めた XVA のボラ感応度は、バンプとは「著しく異なり」、「重大なミスヘッジにつながる」（§4.2、図5）。誤差の大きさは図で示されるだけで、表にはなっていない。
- 回帰の adjoint の計算方法:
  - 前進計算は β=Ψ⁻¹Ω（Ψ は正規方程式の行列、Ω は目的変数に基底を掛けたもの）。
  - 逆方向は Ω̄=Ψ⁻ᵀβ̄、Ψ̄=−Ω̄βᵀ で（式(75)、Giles 2008 の行列微分則による）、得た Ω̄ と Ψ̄ を基底関数、次の時点の価値、ニューメレール、状態変数へ伝播する（式(76)–(78)）。
  - 追加コストは基底次元の小さな線形代数で済む。
- 出典: https://discovery.ucl.ac.uk/id/eprint/1523595/1/Macrina_AADMC_final.pdf
- Capriotti & Giles のレビュー（1.6）も、β 固定は行使境界が最適に近い場合の近似であり、XVA では回帰の寄与を省くと「重大なミスヘッジ」につながり得ると要約している。

**2.5 Fries, "Automatic backward differentiation for American Monte-Carlo algorithms (conditional expectation)", SSRN/arXiv:1707.04942（2017）、*Quantitative Finance* 19(6), 2019, pp.1043–1059** [本文で確認]
- 本文の主な主張:
  - 標準的なバミューダンでは、最適な行使判定のおかげで境界の1次の効果が消えるので、境界を微分する必要はない（§2.2）。これがエンベロープ定理的な議論にあたる。
  - ただしこれはペイオフの形と1次感応度に限った話である。バミューダン・デジタルでは成り立たず、境界を固定するとデルタが0になってしまう（§2.2、§4）。
- 提案手法:
  - 逆方向では「adjoint の条件付き期待値」を取ることで、回帰そのものを微分しない（§3.2、§6 "Do not differentiate an approximation"）。
  - 価値評価と微分で異なる回帰基底を使ってもよい（§6）。
- 出典: https://arxiv.org/pdf/1707.04942 、掲載誌 https://ideas.repec.org/a/taf/quantf/v19y2019i6p1043-1059.html

**2.6 Huge & Savine, "Differential Machine Learning", arXiv:2005.02347（v4、2020年9月）** [確認]
- AAD で得た pathwise 微分をラベルの微分として学習に使う手法。用途として XVA、CCR、FRTB、SIMM-MVA を挙げている。
- 出典: https://arxiv.org/abs/2005.02347
- 将来時点の価格や感応度の代理モデル（MVA 用のSIMM入力など）の候補になる [推測]。

**2.7 その他: Piterbarg と Leclerc 他** [抄録のみ]
- Piterbarg, "Computing deltas of callable LIBOR exotics in forward LIBOR models", *JCF* 7(3), 2004: デルタを最適行使時刻とクーポンのデルタの汎関数として表す。出典: https://www.risk.net/node/2160553
- Leclerc, Liang, Schneider, "Fast Monte Carlo Bermudan Greeks", *Risk* 2009: Giles–Glasserman の手法をバミューダンに拡張。出典: https://www.risk.net/node/1530704
- どちらも行使戦略を固定する系譜とされるが、本文では [未確認]。

---

### 3. Antonov 系と Andreasen

- **Antonov, "Algorithmic differentiation for callable exotics"（Numerix 研究ノート、2016年11月）** [抄録のみ]
  - 回帰がパス間の依存を生むので、パスごとに行う通常の AD は直接適用できないと指摘している。
  - そのうえで、回帰を含めるよう AD を拡張し、テープを使わない backward differentiation（BD）を提案している。
  - 1ファクター Hull-White のバミューダン・スワップションで、50–100個のグリークスを単一価格計算の2–6倍で計算し、バンプより最大20倍速いと報告している。
  - 出典: https://www.numerix.com/node/253
- **Antonov, Issakov, Konikov, McClelland, Mechkov, "Pathwise XVA Greeks for early-exercise products", *Risk* 2018年1月** [ベンダー抄録のみ]
  - 「一部の商品では回帰推定の微分を回避できる」と主張し、価格計算とほぼ同じコストでグリークスを計算できるとしている。回避できる条件は本文未読のため [未確認]。
  - 出典: https://www.numerix.com/resources/quantitative-research/pathwise-xva-greeks-for-early-exercise-products
- **Antonov, McClelland, Issakov, "Efficient SIMM-MVA calculations for callable exotics", *Risk* 2018年9月**（5節参照）。出典: https://www.numerix.com/node/260
- **Andreasen "CVA on an iPad Mini"（2014年頃）** [未確認]
  - 一次資料（スライドや論文）は見つからなかった。
  - Risk.net の "In-house system of the year: Danske Bank"（https://www.risk.net/node/2387869 ）は二次資料にとどまる。
  - Andreasen の xVA 再定式化の内容は、2.3 のセミナー抄録で「未公開」とされている。

---

### 4. 陰関数定理によるキャリブレーション感応度

- **Henrard, "Adjoint algorithmic differentiation: calibration and implicit function theorem", *JCF* 17(4), 2014, pp.37–47（DOI 10.21314/JCF.2014.284）** [抄録のみ]
  - キャリブレーションを含む価格計算の全体を微分する際、IFT を使えばソルバー自体を微分する必要はないと述べている。
  - 出典: https://www.risk.net/journal-of-computational-finance/2347676/
  - 2011年の OpenGamma ワーキングペーパー版は [未確認]。
- **Christianson, "Reverse accumulation and attractive fixed points", *Optimization Methods and Software* 3(4), 1994, pp.311–326** [抄録で確認]
  - 収束する固定点反復を含む関数を逆方向に累積する手法。
  - 勾配の精度を関数値と同等にでき、コストは独立変数の数によらない。反復が逐次でも入れ子でも適用できる。
  - 出典: https://uhra.herts.ac.uk/id/eprint/421/
- **Goloubentsev, Lakshtanov, Piterbarg, "Automatic implicit function theorem", *Risk* 2022年9月** [抄録のみ]
  - exact-fit なキャリブレーションの微分を完全に自動化する手法。
  - 出典: https://www.risk.net/cutting-edge/7954505/automatic-implicit-function-theorem
- 適用上の注意:
  - IFT が成り立つには、キャリブレーションが exact-fit（あるいは最適性条件 ∇=0）を満たし、ヤコビアンが正則であることが前提になる。最小二乗フィットなら最適性条件の2階微分が必要になる。
  - ソルバーが収束しきっていない分は、感応度の誤差になる（Christianson の停止基準の議論を参照）。

---

### 5. MVA/SIMM

- **Green & Kenyon, "MVA: initial margin valuation adjustment by replication and regression", arXiv:1405.0508（2014）、*Risk* 2015 版あり** [本文で確認]
  - MVA は期待 IM プロファイル E_t[I(u)] の積分として表される。
  - VaR 型の IM は LSM 回帰（ポートフォリオ価値の回帰と shocked-state 拡張）で近似している。
  - 本文には SIMM のような感応度ベースの IM は出てこない。ショックは外生的に固定されている。
  - 出典: https://ar5iv.arxiv.org/html/1405.0508
- **Caspers, Giltinan, Lichters, Nowaczyk, "Forecasting initial margin requirements – a model evaluation", *J. Risk Mgmt in Fin. Inst.* 10(4), 2017** [抄録のみ]
  - 回帰型の動的 IM モデルを SIMM と比較し、単一取引では良好だが、分散の効いた大きなポートフォリオでは精度が落ちると報告している。
  - 出典: https://hstalks.com/article/2900/
  - Anfuso 他の動的 IM 論文（"A sound modelling and backtesting framework…"）は、本調査では原典を確かめられていない [未確認]。
- **ISDA SIMM（公開版 v2.9+2606、2026年12月5日発効と記載）** [確認、ただし抽出テキストには崩れあり]
  - リスククラスは6つ、プロダクトクラスは4つ。
  - Delta、Vega、Curvature の各マージンは次の手順で計算する。
    1. 加重感応度 WS = RW × s × CR（CR は集中度ファクター）を作る。
    2. バケット内を ρ で平方根集約する。
    3. バケット間を γ で集約する。
    4. リスククラス間を ψ で集約する。
    5. プロダクトクラス間は単純和を取る。
  - IR デルタの定義は「市場レート（ゼロレートではない）1bp シフトでの PV01」。
  - 出典: https://www.isda.org/a/ZYxiE/ISDA-SIMM_v2.92606_PUBLIC.pdf （v2.7: https://www.isda.org/a/2IOgE/ISDA_SIMM_v2.7_PUBLIC.pdf ）
- **AAD で SIMM-MVA を計算する文献**: Antonov, McClelland, Issakov (2018)（出典は第3節）
  - 将来時点価値の市場クオート感応度を得るため、パラメータ感応度を LSM を通して微分し、それに将来時点のパラメータ→クオートのヤコビアンを掛ける [抄録のみ]。
  - Risk.net の記事（https://www.risk.net/node/5915446 、2018年9月）によれば、シミュレーション部分を AD で、後退評価のパラメータ感応度を tangent で計算し、ブルートフォースより数百倍速い。
  - 依頼文の "Capriotti/Lee AAD for MVA" に当たる原典は見つからなかった [未確認]。

---

### 6. KVA/FVA の非線形性

- **Green & Kenyon, "KVA: Capital Valuation Adjustment", arXiv:1405.0515（v3、2014年10月）** [抄録で確認]
  - Burgard–Kjaer の semi-replication に資本の項を加えたもの。ヘッジ取引自体も規制資本を生むと指摘している。
  - 出典: https://www.arxiv.org/abs/1405.0515
- **Green, Kenyon, Dennis, "KVA: capital valuation adjustment by replication", *Risk* 2014年12月** [本文で確認]
  - 「規制資本はポートフォリオ水準の要件である」としている。
  - 標準的手法の CVA 資本は全カウンターパーティにまたがって計算されると述べている。
  - 新規取引だけの影響を見るのでは不十分で、ヘッジパッケージの増分で価格付けすべきだとしている。
  - 出典: https://discovery.ucl.ac.uk/10140255/1/2014-12%20RISK%20kva.pdf
- ベクトル受け渡しが必要になる理由 [推測。上の記述と SIMM の式構造からの推論]:
  - FVA は、ネッティングセット（あるいはファンディングセット）の純露出の正負で調達・運用スプレッドが変わる非対称な関数 f(Σ_trade V) である。
  - KVA の資本 K(u) は、SA-CCR の EAD、CVA 資本の平方根集約、SIMM の平方根集約などを含み、カウンターパーティやポートフォリオをまたぐ非線形関数になる。
  - そのため、下位ノードがスカラーの XVA を合計しても上位の値にはならない。上位は、パス×時点の露出（あるいは感応度）ベクトルを受け取って非線形関数を評価する必要がある。
  - 逆方向では、上位が ∂F/∂V_{path,t} をシードとして下位に返す。これが合意済みの集約ツリー設計にあたる。
  - Burgard–Kjaer の FVA 原典はこのセッションで再確認していない [未確認]。

---

### 7. PV調整（モデル価格とフロント価格の差）

- 直接の文献は見つからなかった。
- ORE（Open Source Risk Engine）の公式ドキュメントにも、Acadia の公開資料にも、"adj×decay(t)" 型の調整を明示した記述は見当たらなかった [未確認]。
- Deelstra, Grzelak, Wolf, "Accelerated computations of sensitivities for xVA", arXiv:2211.17026 は、市場クオートへのキャリブレーション（Hull-White のドリフトをイールドカーブに合わせる）を論じている。ただし、t=0 の差分を別個に調整する方法は扱っていない。また、カーブ感応度を計算する際に他のモデルパラメータは通常は再キャリブレーションしない、と述べている（付録B）。出典: https://arxiv.org/pdf/2211.17026
- [推測] adj を入力ノードにした場合の adjoint は adj̄ = Σ_t ∂XVA/∂V(t)·decay(t) と単純になる。これにフロント側の ∂adj/∂market（フロントのデルタからモデルのデルタを引いたもの）をつなげるのは連鎖律そのものである。文献上の裏付けがないので、decay の形、定数固定モードとの差、バックテストを社内で検証する必要がある。

---

### 8. bump-over-AAD と AAD-over-AAD（forward-over-reverse）

- **Capriotti & Giles レビュー（2024）** [確認]
  - AAD でヘッセ行列を得るコストは入力数に比例する。
  - AAD 勾配の有限差分でヘッセ行列を推定する方法も同じコストのオーダーで、実装が簡単だという利点がある。
  - 出典: 1.6 に同じ。
- **Naumann & du Toit, "Adjoint algorithmic differentiation tool support for typical numerical patterns in computational finance"（NAG 技術報告、年記載なし。2014年以降と推定）** [確認]
  - ヘッセ行列全体は2次 adjoint コード（tangent-over-adjoint）を n 回実行して得られ、ヘッセ行列×ベクトル積なら1回の実行で済む。
  - forward-over-forward のコストは2次有限差分と同じオーダーで、有限差分の精度は「一般に満足できるものからほど遠い」（特に単精度）とされる。
  - モンテカルロでヘッセ行列を計算した例（n=222）では、中心差分が6.1時間、forward モードが3.2時間かかった。素朴な adjoint はメモリ不足になり、アンサンブル単位でチェックポイントを取った adjoint は13.8分・8.7MB で済んだ。
  - 出典: https://support.nag.com/doc/techrep/pdf/tr3_14.pdf
- **Pagès, Pironneau, Sall, "Vibrato and automatic differentiation for high order derivatives…", arXiv:1606.06143（2016）** [確認]
  - 不連続ペイオフに AD を当てて2次感応度を求めると、Dirac 近似のパラメータに敏感で「危険」だとしている。
  - vibrato と AD の組み合わせを推奨している。
  - 出典: https://arxiv.org/pdf/1606.06143.pdf
- bump-over-AAD と forward-over-reverse を XVA で定量比較した原典は見つからなかった [未確認]。

---

### 目標アーキテクチャへの示唆

1. **回帰と集計をネッティングセット単位の1タスクに統合（順方向・逆方向を同一タスク内）**
   - 裏付け: 回帰の adjoint はパスをまたぐリダクション（Ψ̄、Ω̄）を必要とする（Antonov 2016、Capriotti–Jiang–Macrina 2017）。全パスが同じタスクにあれば自然に実装できるので、統合方針と整合する。
   - リスク: ネッティングセットの大きさが偏ると、タスクの負荷も偏る。
2. **グローバルバリアの廃止と集約ツリー（非線形集約ではベクトルを渡し、逆方向は seed を配る）**
   - 裏付け: KVA の資本はポートフォリオ水準の要件である（Green–Kenyon–Dennis 2014）。SIMM は平方根集約や集中度ファクターを含む非線形関数である（ISDA SIMM）。どちらもスカラーの合算では扱えないので、ベクトル受け渡しが必要になる [推測を含む]。
   - リスク: 上位ノードがパス×時点ベクトルを保持するメモリとネットワーク帯域。CVA 資本がカウンターパーティをまたぐことで、ツリーの根に負荷が集中すること。
3. **逆方向でβを定数として固定**
   - 文献上の妥当性は**限定的で、最大のリスク**になる。
   - エンベロープ定理的な議論が成り立つのは、回帰が最適行使の判定だけに使われ、かつ1次感応度を求める場合に限られる（Fries 2017/2019 §2.2、Capriotti & Giles 2024）。デジタル型では破綻する（Fries）。
   - 回帰した継続価値をそのまま露出（EE や V⁺）として使う場合は、β 固定は「まったく使えない」とされ、ミスヘッジが生じたと報告されている（Capriotti–Jiang–Macrina 2017 §3.4、§4.2）。
   - Savine/Huge の「LSM Reloaded」が β を固定しているかどうかは、原典が非公開のため [未確認]。Antonov 他 (2018) が回帰の微分を回避できるとする条件も [未確認]。
   - 対策案 [推測]:
     - (a) 回帰は行使判定と露出の符号判定（指示関数）だけに使い、露出の値は実現キャッシュフローから組み立てる（Andreasen 型と言われるが [未確認]）。こうすれば準最適性の議論が使える余地がある。
     - (b) β の adjoint（式(75)、正規方程式の逆行列を1回適用するだけで安価）を実装し、固定モードはオプションにする。
     - (c) Fries の「adjoint の条件付き期待値」方式を使う。
     - どれを選ぶ場合でも、バンプとの比較による検証を必須にする。
4. **チェックポイント**
   - 裏付け: Capriotti & Giles (2012)、Savine (2018) 第13章、Naumann & du Toit（チェックポイントでメモリが数百MB→MB級に減った例）。
   - 注意: 回帰の逆方向ではパス全体の adjoint が必要になるので、時点ごとのチェックポイントと相性がよい。
5. **ハイブリッドモデルと、線形商品の解析式**
   - 解析式を使っても AAD とは矛盾しない（Capriotti & Lee 2014 は半解析モデルにも AAD を適用している）。
   - 回帰をエキゾチックに限定すれば、上記3のリスクがかかる範囲も絞れる。
6. **キャリブレーションを IFT で微分**
   - 裏付け: Henrard 2014、Christianson 1994、Capriotti–Jiang–Macrina の AAD-PDE（約3桁の高速化）、Goloubentsev 他 2022。
   - リスク: exact-fit でない場合（最小二乗フィットなど）の扱いと、ソルバーの収束許容誤差。
   - SIMM-MVA には将来時点のパラメータ→クオートのヤコビアンも必要になる（Antonov–McClelland–Issakov 2018）。SIMM の IR デルタは市場レートに対する PV01 と定義されている。
7. **PV調整（adj の入力ノード化と定数固定モード）**
   - 直接の文献はない。連鎖律としては自明だが、規制や監査の観点で裏付けがない設計判断になるので、社内で文書化し検証する必要がある。
8. **クロスガンマを bump-over-AAD で計算**
   - 裏付け: Capriotti & Giles 2024 は、AAD 勾配の有限差分を tangent-over-adjoint と同じコストのオーダーで、実装が簡単な方法として挙げている。
   - リスク: 有限差分の精度（Naumann & du Toit）。不連続点での2次感応度の不安定さ（Pagès 他）は、指示関数の平滑化（Savine & Andreasen 2021 第IV部の fuzzy logic）とセットで扱う必要がある。
   - 日次クロスガンマなら、バンプ数×AAD のコストで足りる。イントラデイ XVA リスクの目標には、Differential ML による代理モデル（Huge & Savine 2020）が補完策の候補になる [推測]。

---

# 第3部 実装調査ノート


凡例: 【ベンダー公表値】は第三者が独立に確認していない数値。【未確認】は一次資料で確認できなかった事項。

---

### 0. 設計に直結する重要な発見（先に要点）

1. **「逆伝播でβ固定」はXVAでは誤差源になりうる。** Capriotti, Jiang, Macrina は、Bermudan では行使境界が準最適なら回帰係数の随伴を無視してよいが、XVAでは「準最適性の議論が使えない」ため正当化できないと述べている。数値例では、β固定AADのボラ感応度がバンプと大きく乖離し、β込みのAADはバンプと一致した。ハザードレート感応度は影響を受けない。全感応度のコストはXVA評価1回の約3倍で、バンプは20倍超（論文の数値）。同論文 §3.4 に回帰ステップの完全な随伴（係数・基底・状態ベクトル）がある。<https://discovery.ucl.ac.uk/id/eprint/1523595/1/Macrina_AADMC_final.pdf>
   → 推奨: βの随伴を正しく扱う（正規方程式／QRの随伴）か、少なくとも検証で「β固定」と「β込み」の差をバンプと比べて定量化すること。
2. **ORE の計算グラフ（CG）では、条件付き期待値の勾配が未実装。** `randomvariable_ops.cpp` の ConditionalExpectation の grad は `QL_FAIL("gradient of conditional expectation not implemented")` になっている。OREがXVA感応度でこれをどう回避しているかは【未確認】。<https://raw.githubusercontent.com/OpenSourceRisk/Engine/master/QuantExt/qle/math/randomvariable_ops.cpp>

---

### 1. C++ AADライブラリ比較

| ライブラリ | 方式 | 前進/ベクトルモード | チェックポイント・外部関数 | スレッド | ライセンス | QuantLib/ORE |
|---|---|---|---|---|---|---|
| **XAD** (Xcelerit) | 演算子オーバーロード＋テープ。任意のJITバックエンドで「1回記録・何度も再生」 | 前進・随伴とも任意階数。ベクトルモードあり | checkpointing と external function interface あり | アクティブテープはスレッドローカル。スレッドごとにテープを持つ | AGPL-3.0。本番利用は有償サブスクリプションが要る場合あり。xad-codegen は商用 | QuantLibAAD / QuantLib-Risks(-Py) あり。ORE連携は【未確認】 |
| **AADC** (MatLogica) | オーバーロード型（`idouble`）で1回トレースし、x64 AVX2/AVX-512 機械語へストリーミングJIT。随伴も同時に生成 | カーネルはベクトル化・マルチスレッド | 分岐は Branch Manager と `iIf()`/`condAssign()` で関数化。ソースのない関数はラッパ経由で局所微分またはローカルバンプ。AIFTでチェックポイント使用 | コンパイル後のカーネルは「元コードが非スレッドセーフでもスレッドセーフ」（ベンダー記載） | 商用（Enterprise／Desk）。Python Community Edition は非商用・学術向けに無償。C++ SDK と AVX-512 は Enterprise | QuantLib・ORE を型置換で統合（ベンダー記載）。GPU出力は無し（NVIDIA倍精度は開発中） |
| **CoDiPack** | 式テンプレート（オーバーロード）。複数のテープ実装: Jacobian/Primal × linear/index | RealForward あり | external functions あり | 【未確認】（関連プロジェクト OpDiLib は未調査） | GPLv3。その他ライセンスは要相談 | 【未確認】 |
| **Adept 2** | 式テンプレート、随伴。配列・行列演算も微分可能 | 【未確認】 | 【未確認】 | 【未確認】 | Apache-2.0（ライブラリ部分） | 無し【未確認】 |
| **dco/c++** (NAG) | ヘッダオンリーのオーバーロード型。v4 でコード生成（オーバーロードとソース変換のハイブリッド）を追加 | ベクトル化タンジェント／随伴。v5.0 で TensorAD | ユーザ定義テープコールバックで手書き・記号的随伴（IFT等）と結合。専用 revolve API は【未確認】 | 共有メモリ・分散メモリ並列をサポート（詳細は【未確認】） | 商用 | Scotiabank のXVAエンジンでAADライブラリを提供（後述） |
| **Enzyme** | LLVM/MLIR IR 上のソース変換 | 逆モード（前進モードは公式トップに明記なし、【未確認】） | 【未確認】 | GPUカーネルの逆モードAD論文あり | Apache-2.0 with LLVM Exceptions | 無し【未確認】 |
| **CppAD** | オーバーロードで ADFun に記録し、再生 | 前進・逆とも | `chkpoint_two`、atomic 関数 | `thread_alloc`/`parallel_setup` 系API、OpenMP等の例あり | EPL-2.0 または GPL-2.0+ | 無し【未確認】 |
| **autodiff** | 前進は `dual`、逆は `var`（式ツリー） | 前進が中心 | 【未確認】 | 【未確認】 | MIT | 無し |

出典:
- XAD: <https://github.com/auto-differentiation/xad>, <https://auto-differentiation.github.io/>（ライセンス、JIT、「XAD比 4〜17倍」「バンプ比 6.5〜83倍」はいずれも【ベンダー公表値】）
- XADテープAPI: <https://auto-differentiation.github.io/ref/tape/>
- QuantLib-Risks の例: 69感応度で評価時間は1.87倍、バンプ比37.4倍【ベンダー公表値】 <https://auto-differentiation.github.io/quantlib-aad/practical-guide>
- AADC: <https://matlogica.com/>, <https://matlogica.com/resources/faq>, <https://www.intel.com/content/www/us/en/developer/articles/technical/accelerate-simulations-aadc-python-cpp-analytics.html>。「6〜1000倍」は【マーケティング文言】。xVA例の「35秒→1秒」は、評価のみの基準と評価＋感応度との比較で同条件ではない（Intel記事が明記）。ORE「100万FXトレードの再計算 98秒→0.4秒（記録350秒）」は【ベンダー公表値】。
- CoDiPack: <https://github.com/SciCompKL/CoDiPack> / Adept: <https://github.com/rjhogan/Adept-2> / dco/c++: <https://www.nag.com/automatic-differentiation>（「最大36,000倍」は【マーケティング文言】）/ Enzyme: <https://github.com/EnzymeAD/Enzyme>, <https://enzyme.mit.edu/>, ライセンスは <https://raw.githubusercontent.com/EnzymeAD/Enzyme/main/LICENSE> / CppAD: <https://github.com/coin-or/CppAD>, <https://cppad.readthedocs.io/latest/multi_thread.html> / autodiff: <https://github.com/autodiff/autodiff>

**XADのAPIで今回の設計に効くもの**（<https://auto-differentiation.github.io/ref/tape/>）
- `clearAll()`: 内容は消えるが確保済みメモリは残る。ドキュメントはパス単位MC-ADでの再利用を想定している。
- `getPosition()`、`resetTo(pos)`、`clearDerivativesAfter(pos)`、`computeAdjointsTo(pos)`: マーク＆巻き戻しに使う。
- `clearDerivatives()`: 記録はそのままで、別の出力のシードで逆伝播し直せる。
- `newNestedRecording()`、`insertCallback`/`pushCallback`、`getAndResetOutputAdjoint`/`incrementAdjoint`: チェックポイント／カスタム随伴。LSMの随伴やネッティング集約ノードを手書きする受け口になる。
- `getMemory()`: テープのバイト数。容量見積りの実測に使う。
- CMakeオプション `XAD_TAPE_REUSE_SLOTS`、`XAD_TAPE_TRACK_LIVE_SLOTS`: 記録を繰り返すときのスロット再利用。

---

### 2. Open Source Risk Engine (ORE)

- **AMC**: v8 でBermudanエクスポージャ用のAMC部品を追加し、v9 でAMCシミュレーションを統合（Example 39）。v13 では Close-out 日の再キャリブレーション、Equity、パス／回帰モデルのシリアライズを追加。v17 では AMC-CG の回帰性能改善、スクリプトトレードのAMC全面対応、AMC学習用の CUDA版 `qrSolve()` を追加。<https://www.github.com/OpenSourceRisk/Engine/blob/master/News.txt>
- **AAD／XVA感応度**: v12 でAADによるXVA感応度のPoC（Example 56）と、Example 61「Fast Sensitivities using AAD and GPUs」を追加。v13 で「AAD based xva sensitivities and dynamic SIMM を支える AMC exposure framework」を追加（同 News.txt）。LSEG公開のロードマップでも、AAD-XVA感応度は「Example 56 で PoC 段階」とされている。<https://www.lseg.com/content/dam/post-trade/en_us/documents/post-trade-solutions/reports/contributions-to-open-source-risk-engine-roadmap.pdf>
- **計算グラフ（CG）と外部計算デバイス**: v11 で外部計算デバイスI/F（GPU）、v12 で OpenCL の倍精度オプションと Mersenne Twister を追加（News.txt）。ソースの構成は以下のとおり。
  - `QuantExt/qle/ad/computationgraph.hpp` に `cg_add`/`cg_mult`/`cg_max`/`cg_conditionalExpectation`/`cg_indicatorGt` などがある。並列化単位とみられる red block もある。<https://raw.githubusercontent.com/OpenSourceRisk/Engine/master/QuantExt/qle/ad/computationgraph.hpp>
  - `backwardderivatives.hpp` は逆伝播を担い、ノード削除（keepNodes）でメモリを抑える。
  - `computeenvironment.hpp` の `ComputeContext` は `createInputVariable`→`applyOperation`→`finalizeCalculation` という流れで、Settings に `useDoublePrecision` と `regressionOrder` がある。<https://raw.githubusercontent.com/OpenSourceRisk/Engine/master/QuantExt/qle/math/computeenvironment.hpp>
  - CMake に `ORE_ENABLE_OPENCL` と `ORE_ENABLE_CUDA` がある（`cmake/commonSettings.cmake`）。
  - `OREAnalytics/orea/engine/xvaenginecg.hpp/.cpp` の `XvaEngineCG` は、CVAノードを「感応度を取る最終結果」として CG 上に構築する。CVAの式は Σ PD_i·max(V_i,0)。オプションに `bumpCvaSensis`、`useExternalComputeDevice`、`useRedBlocks`、`regressionOrder`、dynamic IM がある。<https://raw.githubusercontent.com/OpenSourceRisk/Engine/master/OREAnalytics/orea/engine/xvaenginecg.cpp>
  - 条件付き期待値の勾配は未実装（§0参照）。
- **XAD／AADCとの統合**: ORE本体のビルドオプションに XAD・AADC は見当たらない【未確認】。ORE本体は自前のCG／AD実装を使っている。AADCはORE向けにブラックボックス統合（`CSVLoader::loadFile` などをフック）があると主張している【ベンダー公表】<https://matlogica.com/resources/faq>。
- **SIMM**: v10 で ISDA SIMM（全バージョン）、v11 で 2.6、v12 で CRIFからの 2.7 を追加。v13 では Dynamic SIMM 用の IR/FX CRIF、v17 では stressed SIMM と「SIMM 2.8+2512」を追加（News.txt）。
- 注: <https://opensourcerisk.pages.dev/roadmap> は公式ミラーかどうか【未確認】のため、補助的にのみ使用した。

---

### 3. メモリとテープ設計

- **パス単位テープの再利用（Savine流）**: Savine の公開コード CompFinance の流れは次のとおり。初期化の直後に `tape.mark()` を打ち、各パスで記録→`propagateToMark()`→`rewindToMark()` を繰り返す。最後にマークより前（モデル初期化部）を `propagateMarkToStart()` で1回だけ逆伝播する。マルチスレッド版はスレッドごとのテープで同じことをする。<https://raw.githubusercontent.com/asavine/CompFinance/master/mcBase.h>（書籍の構成は 10.2 "Memory Management and the Tape Class"、13.1 "Check-Pointing" <https://oreilly.com/library/view/modern-computational-finance/9781119539452/c10.xhtml>）。XADでは `getPosition`/`resetTo`/`computeAdjointsTo` で同じ形を取れる。
- **今回の設計への当てはめ**: LSMとFVA/KVAの非線形集約は「全パスの横断情報」を要する。そのため、純粋な1パス完結型にはできない。現実的な形は次の3段になる。(a) 順方向で全パスのパス×時点ベクトルを保持する。(b) 回帰・集約ノードは**ベクトル単位の外部関数（カスタム随伴）**として実装する。(c) パスごとの局所計算は、パス単位で記録→逆伝播→巻き戻し、または再計算で処理する。ORE の CG／RandomVariable 設計（ノードの値がパス数長ベクトル）も同じ発想である。
- **チェックポイント**: 時点チェックポイントでは、時点 t_k の状態（短期金利、FX、LGM状態など）だけを保存し、逆伝播時にその区間を再計算する。最適な配置は Griewank & Walther の binomial checkpointing（Algorithm 799 "revolve", ACM TOMS 2000）が与える。<https://bibbase.org/network/publication/griewank-walther-algorithm799revolveanimplementationofcheckpointingforthereverseoradjointmodeofcomputationaldifferentiation-2000>。JAXのドキュメントも、再帰的チェックポイントで保存量を O(log D)、FLOPを O(log D) 倍にできると説明している。<https://docs.jax.dev/en/latest/gradient-checkpointing.html>
- **保持コストの見積り方**: パス×時点ベクトル1本のバイト数は N_path × N_time × 8B。例えば 1万パス×200時点なら約16MB で、随伴を持つと倍になる。これに、上位へ渡す本数（ネッティングセットごとにV、担保後E、資金調達用など）と同時に生きているノード数を掛ける。テープ側は「1パスあたりの記録バイト数（XADなら `getMemory()` で実測）× 同時記録パス数（≒スレッド数）」で見積もる。dco の例は GPUで約420MB、CPUテープで約268MB（1万パス×360ステップ、438入力）。保存はlog価格とローカルボラに限定し、残りは再計算していた。<https://support.nag.com/doc/techrep/pdf/adjoint-algorithmic-differentiation-of-gpu-accelerated-app.pdf>
- **テープ圧縮**: CoDiPack の index／primal-value テープは Jacobian テープより小さい（README）。dco は局所ヤコビアンの前累積（preaccumulation）でテープを縮小する（NAGページ）。AADCは「コード折り畳みと圧縮」でフットプリントを抑えるとしており、LSMのAADでもメモリ増は約2倍にとどまるという【ベンダー公表値】。
- **AADC の「記録1回・再生多数」**: 1回のトレースから順方向カーネルと随伴カーネルを生成し、パスごとに再生する。確率的分岐は `iIf` などで関数化が必須で、取りこぼした分岐は追跡される（Intel記事、FAQ）。分岐の多いペイオフは書き換えコストが大きい。

---

### 4. GPU／JAX との比較

- **JAX**: `jax.checkpoint`（remat）は保存する残差を制御し、ポリシー（`save_only_these_names` など）やホストへのオフロードを選べる。`lax.scan` の本体に掛けることが推奨されている。<https://docs.jax.dev/en/latest/gradient-checkpointing.html>。**デフォルトはfloat32**で、`jax_enable_x64` を有効にしないと倍精度にならない。<https://docs.jax.dev/en/latest/101/default_dtypes.html>。vmap／XLA の詳細は本調査では精読していない【未確認】。
- **GPUの倍精度**: H100 SXM の FP64 は 34 TFLOPS（Tensor Core で 67）、HBM 80GB、帯域 3.35TB/s【ベンダー公表値】。<https://www.nvidia.com/en-us/data-center/h100.html>。テープとパス×時点ベクトルはGPUメモリ容量に制約される。
- **事例**:
  - NAG・RWTH: CUDA MCカーネルの随伴を手書きし、dco の external function でテープに接続した。全勾配で 522ms、CPU単スレッドのタンジェントでは2時間超（論文値）。<https://support.nag.com/doc/techrep/pdf/adjoint-algorithmic-differentiation-of-gpu-accelerated-app.pdf>
  - Scotiabank: Azure GPU上でAAD（NAGライブラリ）とDNNを組み合わせた。GPUはCPUの10〜50倍で、XVA価格は約20秒（Risk.net記事、【ベンダー/当事者公表値】）。<https://www.scotiabank.com/content/dam/scotiabank/corporate/news/assets/Technology_innovation_of_the_year_Scotiabank_Risknet.pdf>
  - Abbas-Turki, Crépey, Diallo: GPUでのネステッドMCによるCVA/FVA/KVA。<https://math.maths.univ-evry.fr/crepey/papers/XVA-NMC-GPU.pdf>
  - ORE: OpenCL、CUDAの `qrSolve`（AMC学習）（News.txt）。
  - GPUでのLSMのバッチ最小二乗: ORE の CUDA qrSolve が一次資料として確認できた唯一の事例。他の銀行事例は【未確認】。
- **CPU AAD と GPU の得失**
  - 帯域・スループット: GPUは帯域とFP64性能で優位。ただしそれを活かすには、データ並列で分岐の少ないカーネルであることが前提になる。
  - テープ容量: CPUは大容量DRAMを使え、テープを置きやすい。GPUは80〜94GBが上限で、チェックポイント／再計算が前提になる。
  - 分岐の多いペイオフ: CPUのオーバーロード型AADがもっとも素直に扱える。AADCとGPUはどちらも分岐の関数化が必要。
  - C#統合: C++のCPUライブラリはC APIで単純に包める。JAXはPythonランタイムが要り、C#から使うにはプロセス外サービスかXLA AOT（未調査）が必要。
  - 精度: JAXはfloat32が既定。

---

### 5. C API／P/Invoke 設計（Microsoft公式）

- **LibraryImport（.NET 7〜）**: コンパイル時にマーシャリングコードを生成し、実行時のILスタブが要らなくなる。Native AOTやトリミングにも適合する。`static partial` と `AllowUnsafeBlocks` が必要。呼出規約は `UnmanagedCallConvAttribute` で指定し、文字列は `StringMarshalling`（UTF-8 も選べる）で指定する。<https://learn.microsoft.com/en-us/dotnet/standard/native-interop/pinvoke-source-generation>
- **例外**: アンマネージ例外の相互運用は**Windowsのみ**サポートされる。Unix ABIには例外処理の定義がないため、未定義動作やクラッシュの原因になる。<https://learn.microsoft.com/en-us/dotnet/standard/native-interop/exceptions-interoperability>
  → C境界のすべてのエクスポート関数で `try{...}catch(...){}` を行い、エラーコードとスレッドローカルの最終エラーメッセージ取得関数を返す（これは一般的な慣行で、MS文書は暗に支持している程度）。
- **不透明ハンドルと所有権**: `SafeHandle` で寿命を管理し、ファイナライザは避ける（DO/AVOID）。<https://learn.microsoft.com/en-us/dotnet/standard/native-interop/best-practices>。C側は `xva_engine_create`/`_destroy` の対で所有権を明示し、C#側では `SafeHandle.ReleaseHandle` で destroy を呼ぶ。
- **配列**: double[]・int[]・blittable構造体はコピー不要。`bool` は非blittableで、C++の1バイトbool とWindows BOOL（4バイト）の不一致に注意（`int` か `[MarshalAs(U1)]` を使う）。ネイティブ側が呼出後もポインタを保持する場合は、`GCHandle` のPinnedまたは `fixed` を使う（同上）。`fixed` は `Span<T>` の `GetPinnableReference` 経由で使えるが、ピン留めは文の本体の間に限られる。<https://learn.microsoft.com/en-us/dotnet/csharp/language-reference/statements/fixed>
  → パス×時点の大きな結果はネイティブ側で所有し、C#からは `Span<double>` で一時的に読むか、呼出側が確保したバッファへコピーさせる設計が安全。
- **コールバック**: デリゲートは `GC.KeepAlive` するか static で保持する（同上）。
- **ABI安定性・スレッド**: C ABIのみを公開し、構造体には version/size フィールドを持たせる（一般慣行で、MS文書に直接の記載はない）。XADのテープはスレッドローカルなので、C#のスレッドプールからではなく、ネイティブ側の自前スレッドプールでタスク（ネッティングセット）を実行する方が安全（XADドキュメントからの推論）。

---

### 6. 検証方法

- **Taylor test**: |J(m+hδm)−J(m)−h∇J·δm| が O(h²) で減ることを確認し、h を半分にして収束率≈2 を見る。<https://dolfin-adjoint-doc.readthedocs.io/en/latest/documentation/verification.html>
- **dot-product test**（タンジェントと随伴の整合 ⟨ẏ, ȳ⟩ = ⟨ẋ, x̄⟩）: 一般手法として採用を推奨する。今回は該当文書を精読していない【出典未確認】。tlm-adjoint の verification モジュールが参考になりうる。<https://tlm-adjoint.github.io/autoapi/tlm_adjoint/verification/index.html>
- **バンプ比較**: 同一乱数シード（共通乱数）で中心差分を取る。ステップは h を対数掃引し、打切り誤差と丸め誤差の谷を選ぶ。max(·,0) や指示関数は不連続でバンプとAADが乖離しやすいので、スムージング（ORE の CG は `eps` を indicator に使用、`xvaenginecg.cpp`）も比較対象にする。
- **β固定と完全随伴の差の定量化**: Capriotti et al. の手順を再現する（§0）。
- **ライブラリの検証ユーティリティ**: XAD・dco/c++ に専用の勾配チェックユーティリティがあるかは【未確認】。dco は v4.1 でテープの CSV/DOT 出力、AADC は `AADC_TRACE=1` トレーサを持つ（NAG/MatLogicaページ）。

---

### 推奨と判断材料

**XVAコアの候補順位**
1. **XAD**: OSSで中身を検証でき、テープAPIが本設計に合う（`clearAll`/`resetTo`/`computeAdjointsTo`/nested recording/callback）。スレッドローカルのテープもタスク並列と相性が良い。QuantLib統合の実績もある。注意点は、AGPLのため社内利用でも商用ライセンスの法務確認が要ることと、性能値がベンダー公表値であること。
2. **AADC**: 記録1回・再生多数でSIMD化とマルチスレッド化が自動で効き、イントラデイ再計算（シナリオ・取引追加の再記録不要）に強い。ただし商用で、C++ SDKはEnterpriseのみ。分岐を `iIf` に書き換える必要があり、GPUにも未対応。XADで設計を固めたあと、ボトルネックのカーネルに試す二段構えが妥当。
3. **dco/c++**: 商用サポートがあり、外部随伴やGPU連携の文献も厚い。銀行XVAでの実績（Scotiabank）もある。費用と、コード開示の制約が判断材料になる。
4. **CoDiPack**（GPLv3）、**CppAD**（EPL/GPL、チェックポイントとマルチスレッドの基盤あり）: 研究・PoC向け。
5. **Enzyme**: ソース変換で高速になりうるが、LLVMツールチェーンへの依存とC#配布の複雑さがある。**Adept**／**autodiff**: XVA規模の実績は【未確認】。

どれを選んでも、**LSM・集約ノードはベクトル単位のカスタム随伴として自作**することになる。ライブラリ差は主に、パスごとの局所テープの速度とメモリに現れる。

**GPU/JAXを採用すべき条件**
- ペイオフとモデルが分岐の少ないデータ並列形で、パス×時点の全状態がGPUメモリ（80〜94GB級）に、チェックポイント込みで収まること。
- FP64性能のあるデータセンターGPU（H100級）を確保でき、JAXでは `jax_enable_x64` を前提にすること。
- C#からの呼出をプロセス外サービス（gRPC等）で許容できること。もしくは CUDA+手書き随伴を C API で包む工数を確保できること。
- 日次クロスガンマ（バンプ×AAD）で「シナリオ数×フルAAD」がCPUグリッドで時間内に終わらない規模であること。これに当たらなければ、CPU AAD（XAD→必要ならAADC）が保守性とC#統合の点で優位。
