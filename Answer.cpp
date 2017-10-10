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

struct targetting_t {
    enum { NONE = -1, DELIVERED = -2 };
    array<int, Parameter::UFOCount> ufo_to_house;
    array<int, Parameter::MaxHouseCount> house_to_ufo;
    int from_ufo(int ufo_index) {
        return ufo_to_house[ufo_index];
    }
    int from_house(int house_index) {
        return house_to_ufo[house_index];
    }
    bool is_targetting(int ufo_index) {
        return ufo_to_house[ufo_index] != NONE;
    }
    bool is_delivered(int house_index) {
        return house_to_ufo[house_index] == DELIVERED;
    }
    void link(int ufo_index, int house_index) {
        unlink_ufo(ufo_index);
        unlink_house(house_index);
        ufo_to_house[ufo_index] = house_index;
        house_to_ufo[house_index] = ufo_index;
    }
    void unlink_ufo(int ufo_index) {
        int & house_index = ufo_to_house[ufo_index];
        if (house_index == NONE) return;
        house_to_ufo[house_index] = NONE;
        house_index = NONE;
    }
    void unlink_house(int house_index) {
        int & ufo_index = house_to_ufo[house_index];
        if (ufo_index == NONE or ufo_index == DELIVERED) return;
        ufo_to_house[ufo_index] = NONE;
        ufo_index = NONE;
    }
    void deliver_house(int house_index) {
        int & ufo_index = house_to_ufo[house_index];
        ufo_to_house[ufo_index] = NONE;
        ufo_index = DELIVERED;
    }
    void clear() {
        fill(whole(ufo_to_house), NONE);
        fill(whole(house_to_ufo), NONE);
    }
    void debug(Stage const & stage) {
        repeat (ufo_index, Parameter::UFOCount) {
            cerr << ufo_to_house[ufo_index] << ' ';
        }
        cerr << endl;
        int house_count = stage.houses().count();
        repeat (house_index, house_count) {
            cerr << house_to_ufo[house_index] << ' ';
        }
        cerr << endl;
    }
};
targetting_t target;

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
    target.clear();
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
    int house_count = stage.houses().count();
    array<int, Parameter::UFOCount> item_count;
    repeat (ufo_index, Parameter::UFOCount) {
        auto const & ufo = stage.ufos()[ufo_index];
        item_count[ufo_index] = ufo.itemCount();

        if (Util::IsIntersect(ufo, stage.office())) {
            actions.add(Action::PickUp(ufo_index));
        }

        if (item_count[ufo_index] != 0) {
            int nearest_house_index = -1;
            double nearest_house_distance = INFINITY;
            repeat (house_index, house_count) if (not target.is_delivered(house_index)) {
                auto const & house = stage.houses()[house_index];

                if (target.from_ufo(ufo_index) == house_index and Util::IsIntersect(ufo, house)) {
                    item_count[ufo_index] -= 1;
                    actions.add(Action::Deliver(ufo_index, house_index));
                    target.deliver_house(house_index);
                    continue;
                }

                if (target.from_house(house_index) == targetting_t::NONE) {
                    double dist = ufo.pos().dist(house.pos());
                    if (dist < nearest_house_distance) {
                        nearest_house_distance = dist;
                        nearest_house_index = house_index;
                    }
                }
            }

            if (not target.is_targetting(ufo_index)) {
                if (item_count[ufo_index] and nearest_house_index != -1) {
                    target.link(ufo_index, nearest_house_index);
                }
            }

            if (item_count[ufo_index] and target.from_ufo(ufo_index) == targetting_t::NONE) {
                repeat (other_ufo_index, Parameter::UFOCount) if (target.is_targetting(other_ufo_index)) {
                    auto const & other_ufo = stage.ufos()[other_ufo_index];
                    int house_index = target.from_ufo(other_ufo_index);
                    auto const & house = stage.houses()[house_index];
                    double this_time = ufo.pos().dist(house.pos()) / ufo.maxSpeed();
                    double other_time = other_ufo.pos().dist(house.pos()) / other_ufo.maxSpeed();
                    if (this_time < other_time) {
                        target.unlink_ufo(other_ufo_index);
                        target.link(ufo_index, house_index);
                        break;
                    }
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
            int house_index = target.from_ufo(ufo_index);
            if (house_index != targetting_t::NONE) {
                target_positions.add(stage.houses()[house_index].pos());
            } else {
                target_positions.add(ufo.pos());
            }
        }
    }
// target.debug(stage);
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
