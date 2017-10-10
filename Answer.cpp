//------------------------------------------------------------------------------
/// @file
/// @author   ハル研究所プログラミングコンテスト実行委員会
///
/// @copyright  Copyright (c) 2017 HAL Laboratory, Inc.
/// @attention  このファイルの利用は、同梱のREADMEにある
///             利用条件に従ってください。
//------------------------------------------------------------------------------

#include "Answer.hpp"

/// プロコン問題環境を表します。
namespace hpc {

static const int TARGET_HOUSE_NONE = -1;

/// 回答コードのサンプルです。
class Solver
{
public:
    void init(const Stage& aStage)
    {
        mTargetHouseIndices.clear();

        // 各UFOは、(houseIndex % ufoCount == ufoIndex)の家を担当する
        int ufoCount = aStage.ufos().count();
        for (int ufoIndex = 0; ufoIndex < ufoCount; ++ufoIndex) {
            mTargetHouseIndices.add(ufoIndex);
        }
    }

    void moveItems(const Stage& aStage, Actions& aActions)
    {
        int ufoCount = aStage.ufos().count();

        for (int ufoIndex = 0; ufoIndex < ufoCount; ++ufoIndex) {
            const auto& ufo = aStage.ufos()[ufoIndex];

            // 農場の上にいたら箱を積み込む
            if (Util::IsIntersect(ufo, aStage.office())) {
                aActions.add(Action::PickUp(ufoIndex));
            }

            // 担当の家の上にいたら配達する
            int houseIndex = mTargetHouseIndices[ufoIndex];
            if (houseIndex == TARGET_HOUSE_NONE) {
                continue;
            }

            if (ufo.itemCount() == 0 || !Util::IsIntersect(ufo, aStage.houses()[houseIndex])) {
                continue;
            }

            aActions.add(Action::Deliver(ufoIndex, houseIndex));

            // 目標の家を更新
            mTargetHouseIndices[ufoIndex] += ufoCount;

            if (mTargetHouseIndices[ufoIndex] >= aStage.houses().count()) {
                // 担当の家はすべて配達終了
                mTargetHouseIndices[ufoIndex] = TARGET_HOUSE_NONE;
            }
        }
    }

    void moveUFOs(const Stage& aStage, TargetPositions& aTargetPositions)
    {
        int ufoCount = aStage.ufos().count();

        for (int ufoIndex = 0; ufoIndex < ufoCount; ++ufoIndex) {
            const auto& ufo = aStage.ufos()[ufoIndex];

            if (ufo.itemCount() == 0) {
                // 箱がなければ農場に向かう
                aTargetPositions.add(aStage.office().pos());
            } else {
                // 箱を積んでいれば担当の家に向かう
                int houseIndex = mTargetHouseIndices[ufoIndex];
                if (houseIndex != TARGET_HOUSE_NONE) {
                    aTargetPositions.add(aStage.houses()[houseIndex].pos());
                } else {
                    aTargetPositions.add(ufo.pos());
                }
            }
        }
    }

private:
    /// 各UFOが次に向かう家のインデックス。
    Array<int, Parameter::UFOCount> mTargetHouseIndices;
};

Solver g_Solver;

//------------------------------------------------------------------------------
/// Answer クラスのコンストラクタです。
///
/// ここに最初のステージの開始前に行う処理を書くことができます。何も書かなくても構いません。
Answer::Answer()
{
}

//------------------------------------------------------------------------------
/// Answer クラスのデストラクタです。
///
/// ここに最後のステージの終了後に行う処理を書くことができます。何も書かなくても構いません。
Answer::~Answer()
{
}

//------------------------------------------------------------------------------
/// 各ステージ開始時に呼び出されます。
///
/// ここで、各ステージに対して初期処理を行うことができます。
///
/// @param[in] aStage 現在のステージ。
void Answer::init(const Stage& aStage)
{
    g_Solver.init(aStage);
}

//------------------------------------------------------------------------------
/// 各ターンの、受け渡しフェーズの行動を指定してください。
///
/// - Actionクラスのstatic関数で行動を生成し、aActionsに格納してください。
/// - 行動は、配列のインデックスの昇順で実行されます。
/// - 同じUFOが複数回行動することもできます。
/// - UFOが箱を持っていないなど、必要な条件を満たしていない場合は何も起こりません。
///   条件の詳細はActionクラスを参照してください。
/// - 1回のフェーズの最大行動数は、Parameter::MaxActionPerTurnです。
/// - aActionsの要素の、UFOや家のインデックスが範囲外の場合、アサートに失敗します。
///
/// @param[in] aStage 現在のステージ。
/// @param[out] aActions この受け渡しフェーズの行動を指定する配列。
void Answer::moveItems(const Stage& aStage, Actions& aActions)
{
    g_Solver.moveItems(aStage, aActions);
}

//------------------------------------------------------------------------------
/// 各ターンの、移動フェーズの行動を指定してください。
/// 
/// - 各UFOごとに、目標座標を指定してください。
/// - aTargetPositions[i]が、aStage.ufos()[i]の目標座標となります。
/// - aTargetPositions.count() == aStage.ufos().count()となるように、配列に座標を追加してください。
///   - 要素数が異なる場合はアサートに失敗します。
/// - aTargetPositionsの要素の値が NaN または ±Inf である場合、アサートに失敗します。
///
/// @param[in] aStage 現在のステージ。
/// @param[out] aTargetPositions 各UFOの目標座標を指定する配列。
void Answer::moveUFOs(const Stage& aStage, TargetPositions& aTargetPositions)
{
    g_Solver.moveUFOs(aStage, aTargetPositions);
}

//------------------------------------------------------------------------------
/// 各ステージ終了時に呼び出されます。
///
/// ここにステージ終了時の処理を書くことができます。何も書かなくても構いません。
///
/// @param[in] aStage 現在のステージ。
void Answer::finalize(const Stage& aStage)
{
}

} // namespace
// EOF
