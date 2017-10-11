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
#include <cassert>
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

struct TargetManager {
    TargetManager() {
        fill(whole(ufo_to_house), NONE);
        fill(whole(house_to_ufo), NONE);
    }
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

int turns_to_move(Vector2 const & from, Vector2 const & to, int max_speed) {
    return ceil(from.dist(to) / max_speed);
}

struct beam_state_t {
    bitset<Parameter::MaxHouseCount> delivered;
    int house;
    int turn;
    vector<int> path;
};
bool operator < (beam_state_t const & s, beam_state_t const & t) { return s.turn < t.turn; }  // weak
bool operator > (beam_state_t const & s, beam_state_t const & t) { return s.turn > t.turn; }
beam_state_t beamsearch(Stage const & stage, bitset<Parameter::MaxHouseCount> const & delivered, UFO const & ufo, int turn_limit) {
    vector<beam_state_t> beam; {
        beam_state_t initial = {};
        initial.delivered = delivered;
        initial.turn = 0;
        initial.house = -1;  // the office
        beam.push_back(initial);
    }
    int house_count = stage.houses().count();
    constexpr int beam_width = 100;
    while (true) {
        vector<beam_state_t> prev_beam;
        prev_beam.swap(beam);
        for (auto const & s : prev_beam) {
            Vector2 pos = s.house == -1 ? stage.office().pos() : stage.houses()[s.house].pos();
            repeat (house_index, house_count) if (not s.delivered[house_index]) {
                beam_state_t t = s;
                t.delivered[house_index] = true;
                t.house = house_index;
                t.turn += turns_to_move(pos, stage.houses()[house_index].pos(), ufo.maxSpeed());
                t.path.push_back(house_index);
                beam.push_back(t);
            }
        }
        sort(whole(beam));
        beam.resize(min<int>(beam.size(), beam_width));
        if (int(beam.front().delivered.count()) == house_count) break;
        if (beam.front().turn > turn_limit) break;
    }
    return beam.front();
}

void move_items_with_large_ufos_plan(Stage const & stage, Actions & actions, TargetManager & target, vector<vector<int> > const & large_path, bitset<Parameter::MaxHouseCount> const & delivered_by_large) {
    int house_count = stage.houses().count();
    array<int, Parameter::UFOCount> item_count;
    repeat (ufo_index, Parameter::UFOCount) {
        auto const & ufo = stage.ufos()[ufo_index];
        item_count[ufo_index] = ufo.itemCount();

        if (item_count[ufo_index] < ufo.capacity() and Util::IsIntersect(ufo, stage.office())) {
            actions.add(Action::PickUp(ufo_index));
        }

        if (item_count[ufo_index] != 0) {
            if (target.is_targetting(ufo_index)) {
                int house_index = target.from_ufo(ufo_index);
                auto const & house = stage.houses()[house_index];
                if (Util::IsIntersect(ufo, house)) {
                    item_count[ufo_index] -= 1;
                    actions.add(Action::Deliver(ufo_index, house_index));
                    target.deliver_house(house_index);
                }
            }

            if (item_count[ufo_index] and not target.is_targetting(ufo_index)) {
                if (ufo.type() == UFOType_Large) {
                    for (int house_index : large_path[ufo_index]) {
                        if (target.from_house(house_index) == TargetManager::NONE) {
                            target.link(ufo_index, house_index);
                            break;
                        }
                    }
                }
            }

            if (item_count[ufo_index] and not target.is_targetting(ufo_index)) {
                int nearest_house_index = -1;
                double nearest_house_distance = INFINITY;
                repeat (house_index, house_count) if (not target.is_delivered(house_index)) {
                    auto const & house = stage.houses()[house_index];
                    if (not delivered_by_large[house_index]) {
                        if (target.from_house(house_index) == TargetManager::NONE) {
                            double dist = ufo.pos().dist(house.pos());
                            if (dist < nearest_house_distance) {
                                nearest_house_distance = dist;
                                nearest_house_index = house_index;
                            }
                        }
                    }
                }
                if (nearest_house_index != -1) {
                    target.link(ufo_index, nearest_house_index);
                }
            }

            if (item_count[ufo_index] and not target.is_targetting(ufo_index)) {
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

void move_ufos(Stage const & stage, TargetPositions & target_positions, TargetManager & target) {
    repeat (ufo_index, Parameter::UFOCount) {
        auto const & ufo = stage.ufos()[ufo_index];

        if (ufo.itemCount() == 0) {
            target_positions.add(stage.office().pos());

        } else {
            int house_index = target.from_ufo(ufo_index);
            if (house_index != TargetManager::NONE) {
                target_positions.add(stage.houses()[house_index].pos());
            } else {
                target_positions.add(ufo.pos());
            }
        }
    }
}

struct turn_output_t {
    Actions actions;
    TargetPositions target_positions;
};
vector<turn_output_t> result;

//------------------------------------------------------------------------------
/// 各ステージ開始時に呼び出されます。
///
/// ここで、各ステージに対して初期処理を行うことができます。
///
/// @param[in] stage 現在のステージ。
void Answer::init(Stage const & a_stage) {
    result.clear();

    for (int turn_limit : { 50, 80, 110 }) {
        Stage stage = a_stage;
        TargetManager target = {};

        vector<vector<int> > path(Parameter::LargeUFOCount);
        bitset<Parameter::MaxHouseCount> delivered_by_large = {};
        repeat (ufo_index, Parameter::LargeUFOCount) {
            auto const & ufo = stage.ufos()[ufo_index];
            assert (ufo.type() == UFOType_Large);
            auto s = beamsearch(stage, delivered_by_large, ufo, turn_limit);
            path[ufo_index] = s.path;
            delivered_by_large = s.delivered;

#ifdef LOCAL
cerr << "path " << ufo_index << ": turn " << s.turn << ": ";
for (int house_index : path[ufo_index]) cerr << house_index << ' ';
cerr << endl;
#endif
        }

        vector<turn_output_t> outputs;
        while (not stage.hasFinished() and stage.turn() < Parameter::GameTurnLimit) {
            turn_output_t output = {};
            move_items_with_large_ufos_plan(stage, output.actions, target, path, delivered_by_large);
            stage.moveItems(output.actions);
            move_ufos(stage, output.target_positions, target);
            stage.moveUFOs(output.target_positions);
            stage.advanceTurn();
            outputs.push_back(output);

#ifdef LOCAL
        // debug
cerr << "turn " << stage.turn() << ": ";
repeat (house_index, stage.houses().count()) cerr << stage.houses()[house_index].delivered();
cerr << " / ";
repeat (ufo_index, Parameter::UFOCount) cerr << target.from_ufo(ufo_index) << "(" << stage.ufos()[ufo_index].itemCount() << ") ";
cerr << endl;
#endif

#ifdef LOCAL
        // check invariant
        repeat (ufo_index, Parameter::UFOCount) {
            int house_index = target.from_ufo(ufo_index);
            if (house_index == TargetManager::NONE) {
                // nop
            } else if (house_index == TargetManager::DELIVERED) {
                assert (false);
            } else {
                assert (target.from_house(house_index) == ufo_index);
            }
        }
        repeat (house_index, stage.houses().count()) {
            int ufo_index = target.from_house(house_index);
            if (ufo_index == TargetManager::NONE) {
                // nop
            } else if (ufo_index == TargetManager::DELIVERED) {
                assert (stage.houses()[house_index].delivered());
            } else {
                assert (target.from_ufo(ufo_index) == house_index);
            }
        }
#endif
        }

        if (result.empty() or outputs.size() < result.size()) {
            result = outputs;
        }
    }
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
    actions = result[stage.turn()].actions;
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
    target_positions = result[stage.turn()].target_positions;
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
