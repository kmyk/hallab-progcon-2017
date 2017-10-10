//------------------------------------------------------------------------------
/// @file
/// @author   ハル研究所プログラミングコンテスト実行委員会
///
/// @copyright  Copyright (c) 2017 HAL Laboratory, Inc.
/// @attention  このファイルの利用は、同梱のREADMEにある
///             利用条件に従ってください。
//------------------------------------------------------------------------------

#include "Answer.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>
#define repeat(i, n) for (int i = 0; (i) < int(n); ++(i))
#define repeat_from(i, m, n) for (int i = (m); (i) < int(n); ++(i))
#define repeat_reverse(i, n) for (int i = (n)-1; (i) >= 0; --(i))
#define repeat_from_reverse(i, m, n) for (int i = (n)-1; (i) >= int(m); --(i))
#define whole(x) begin(x), end(x)
#define unittest_name_helper(counter) unittest_ ## counter
#define unittest_name(counter) unittest_name_helper(counter)
#define unittest __attribute__((constructor)) void unittest_name(__COUNTER__) ()
using ll = long long;
using namespace std;
template <class T> inline void setmax(T & a, T const & b) { a = max(a, b); }
template <class T> inline void setmin(T & a, T const & b) { a = min(a, b); }

/// プロコン問題環境を表します。
namespace hpc {

static const int TARGET_NONE = -1;
static const int TARGET_DELIVERED = -2;
array<char, Parameter::UFOCount> target_house;
array<int, Parameter::MaxHouseCount> targetted_by;

//------------------------------------------------------------------------------
/// Answer クラスのコンストラクタです。
///
/// ここに最初のステージの開始前に行う処理を書くことができます。何も書かなくても構いません。
Answer::Answer() {
}

//------------------------------------------------------------------------------
/// Answer クラスのデストラクタです。
///
/// ここに最後のステージの終了後に行う処理を書くことができます。何も書かなくても構いません。
Answer::~Answer() {
}

//------------------------------------------------------------------------------
/// 各ステージ開始時に呼び出されます。
///
/// ここで、各ステージに対して初期処理を行うことができます。
///
/// @param[in] stage 現在のステージ。
void Answer::init(Stage const & stage) {
    fill(whole(target_house), TARGET_NONE);
    fill(whole(targetted_by), TARGET_NONE);
}

//------------------------------------------------------------------------------
/// 各ターンの、受け渡しフェーズの行動を指定してください。
///
/// - Actionクラスのstatic関数で行動を生成し、actionsに格納してください。
/// - 行動は、配列のインデックスの昇順で実行されます。
/// - 同じUFOが複数回行動することもできます。
/// - UFOが箱を持っていないなど、必要な条件を満たしていない場合は何も起こりません。
///   条件の詳細はActionクラスを参照してください。
/// - 1回のフェーズの最大行動数は、Parameter::MaxActionPerTurnです。
/// - actionsの要素の、UFOや家のインデックスが範囲外の場合、アサートに失敗します。
///
/// @param[in] stage 現在のステージ。
/// @param[out] actions この受け渡しフェーズの行動を指定する配列。
void Answer::moveItems(Stage const & stage, Actions & actions) {
    int ufo_count = stage.ufos().count();
    int house_count = stage.houses().count();
    repeat (ufo_index, ufo_count) {
        auto const & ufo = stage.ufos()[ufo_index];

        if (Util::IsIntersect(ufo, stage.office())) {
            actions.add(Action::PickUp(ufo_index));
        }

        if (ufo.itemCount() != 0) {
            bool is_delivering = false;
            int nearest_house_index = -1;
            double nearest_house_distance = INFINITY;
            repeat (house_index, house_count) {
                auto const & house = stage.houses()[house_index];
                if (house.delivered()) continue;

                if (target_house[ufo_index] == house_index and Util::IsIntersect(ufo, house)) {
                    is_delivering = true;
                    actions.add(Action::Deliver(ufo_index, house_index));
                    target_house[ufo_index] = TARGET_NONE;
                    targetted_by[house_index] = TARGET_DELIVERED;
                    break;
                }

                if (targetted_by[house_index] == TARGET_NONE) {
                    double dist = ufo.pos().dist(house.pos());
                    if (dist < nearest_house_distance) {
                        nearest_house_distance = dist;
                        nearest_house_index = house_index;
                    }
                }
            }

            if (target_house[ufo_index] == TARGET_NONE) {
                if (ufo.itemCount() - int(is_delivering) != 0 and nearest_house_index != -1) {
                    target_house[ufo_index] = nearest_house_index;
                    targetted_by[nearest_house_index] = ufo_index;
                }
            }
        }
    }
}

//------------------------------------------------------------------------------
/// 各ターンの、移動フェーズの行動を指定してください。
/// 
/// - 各UFOごとに、目標座標を指定してください。
/// - target_positions[i]が、stage.ufos()[i]の目標座標となります。
/// - target_positions.count() == stage.ufos().count()となるように、配列に座標を追加してください。
///   - 要素数が異なる場合はアサートに失敗します。
/// - target_positionsの要素の値が NaN または ±Inf である場合、アサートに失敗します。
///
/// @param[in] stage 現在のステージ。
/// @param[out] target_positions 各UFOの目標座標を指定する配列。
void Answer::moveUFOs(Stage const & stage, TargetPositions & target_positions) {
    int ufo_count = stage.ufos().count();
    repeat (ufo_index, ufo_count) {
        auto const & ufo = stage.ufos()[ufo_index];

        if (ufo.itemCount() == 0) {
            target_positions.add(stage.office().pos());

        } else {
            int house_index = target_house[ufo_index];
            if (house_index != TARGET_NONE) {
                target_positions.add(stage.houses()[house_index].pos());
            } else {
                target_positions.add(ufo.pos());
            }
        }
    }
}

//------------------------------------------------------------------------------
/// 各ステージ終了時に呼び出されます。
///
/// ここにステージ終了時の処理を書くことができます。何も書かなくても構いません。
///
/// @param[in] stage 現在のステージ。
void Answer::finalize(Stage const & stage) {
}

} // namespace
// EOF
