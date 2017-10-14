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
#include <memory>
#include <set>
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

/// こちらで勝手に定義した部分はここに
namespace Solver {

/// UFOと家の双方向の対応を管理します。NONE,DELIVEREDなものを除けば全単射が保証されます。
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

}

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

namespace Solver {

/// 座標がstage上にあるかを判定します。
bool is_on_stage(int y, int x) {
    return 0 <= y and y < Parameter::StageHeight and 0 <= x and x < Parameter::StageWidth;
}
bool is_on_stage(Vector2 pos) {
    return is_on_stage(pos.y, pos.x);
}

/// hpc::Stage が実装内部に持ってる定数群です。
namespace StageParameter {
const int MinTownCount = 2;
const int MaxTownCount = 3;
const int TownRadius = 100;
const int TownHouseCount = 20;
const int MinRandomHouseCount = 5;
const int MaxRandomHouseCount = 40;
const int OfficeMargin = 150;
const int FieldMargin = 20;
}

/// 検出された街を表現する構造体
struct town_t {
    Vector2 center;
    vector<int> house_indices;
};
/// 街を検出し列挙します。
vector<town_t> detect_towns(Houses const & houses) {

    // 各位置から半径TownRadiusで見える家の数を数える。imos法でO(Nr + HW)
    typedef array<array<int8_t, Parameter::StageWidth + 1>, Parameter::StageHeight> cnt_array_t;
    auto cnt_ptr = unique_ptr<cnt_array_t>(new cnt_array_t {});
    auto & cnt = *cnt_ptr;
    repeat (house_index, houses.count()) {
        auto const & house = houses[house_index];
        repeat_from (dy, - StageParameter::TownRadius, StageParameter::TownRadius + 1) {
            int y = house.pos().y + dy;
            if (0 <= y and y < Parameter::StageHeight) {
                int dx = ceil(sqrt(pow(StageParameter::TownRadius, 2) - pow(dy, 2)));
                int xl = max(0, min<int>(Parameter::StageWidth - 1, house.pos().x - dx));
                int xr = max(0, min<int>(Parameter::StageWidth - 1, house.pos().x + dx));
                cnt[y][xl] += 1;
                cnt[y][xr + 1] -= 1;
            }
        }
    }
    repeat (y, Parameter::StageHeight) {
        repeat (x, Parameter::StageWidth) {
            cnt[y][x + 1] += cnt[y][x];
        }
    }
#ifdef DEBUG
repeat (y, Parameter::StageHeight) if (y % 10 == 0) {
    repeat (x, Parameter::StageWidth) if (x % 10 == 0) {
        fprintf(stderr, "%d", cnt[y][x]/4);
    }
    fprintf(stderr, "\n");
}
fprintf(stderr, "\n");
#endif

    // 適当な位置からDFSして中心を判断
    int dfs_max_cnt = -1;
    Vector2 dfs_pos;
    function<void (int, int)> dfs = [&](int y, int x) {
        if (dfs_max_cnt < cnt[y][x]) {  // TODO: これだと真の中心にはならない
            dfs_max_cnt = cnt[y][x];
            dfs_pos = Vector2(x, y);
        }
        cnt[y][x] = -1;  // 使用済みflagを同居 汚ないが時空間効率のため
        repeat_from (ny, y - 1, y + 2) {
            repeat_from (nx, x - 1, x + 2) {
                if (is_on_stage(ny, nx) and cnt[ny][nx] >= dfs_max_cnt - 1) {
                    dfs(ny, nx);
                }
            }
        }
    };
    vector<Vector2> town_centers;
    repeat (y, Parameter::StageHeight) {
        repeat (x, Parameter::StageWidth) {
            if (cnt[y][x] >= StageParameter::TownHouseCount) {
                dfs_max_cnt = -1;
                dfs(y, x);
                town_centers.push_back(dfs_pos);
            }
        }
    }

    // 街を復元
    vector<town_t> towns;
    for (auto town_center : town_centers) {
        town_t town = {};
        town.center = town_center;
        repeat (house_index, houses.count()) {
            auto const & house = houses[house_index];
            if (house.pos().dist(town.center) <= StageParameter::TownRadius + 3) { // 3 は余裕
                town.house_indices.push_back(house_index);
            }
        }
#ifdef DEBUG
fprintf(stderr, "town (%d, %d) : size %d : ", int(town.center.y), int(town.center.x), int(town.house_indices.size()));
for (int house_index : town.house_indices) fprintf(stderr, "%d ", house_index);
fprintf(stderr, "\n");
#endif
        towns.push_back(town);
    }

    // 衝突してたら併合
    repeat (town_index, towns.size()) {
        auto const & town = towns[town_index];
        repeat (other_town_index, towns.size()) if (town_index != other_town_index) {
            auto const & other_town = towns[other_town_index];
            if (town.house_indices.size() <= other_town.house_indices.size()) {
                vector<int> intersection;
                set_intersection(whole(town.house_indices), whole(other_town.house_indices), back_inserter(intersection));
                if (intersection.size() >= 10) {
                    towns.erase(towns.begin() + town_index);
                    -- town_index;
                    break;
                }
            }
        }
    }
#ifdef DEBUG
fprintf(stderr, "merged:\n");
repeat (town_index, towns.size()) {
auto const & town = towns[town_index];
fprintf(stderr, "town (%d, %d) : size %d : ", int(town.center.y), int(town.center.x), int(town.house_indices.size()));
for (int house_index : town.house_indices) fprintf(stderr, "%d ", house_index);
fprintf(stderr, "\n");
}
#endif

    // 確認
#ifdef LOCAL
    assert (StageParameter::MinTownCount <= towns.size() and towns.size() <= StageParameter::MaxTownCount);
#endif
    return towns;
}

/// 2点間の移動に要するターン数を計算します。
///
/// @note 改善余地あり。半径とかをまったく考慮していない。
int turns_to_move(Vector2 const & from, Vector2 const & to, int max_speed) {
    return ceil(from.dist(to) / max_speed);
}

/// ビームサーチにおける状態。
struct beam_state_t {
    bitset<Parameter::MaxHouseCount> delivered;
    int house;
    int turn;
    vector<int> path;
};
bool operator < (beam_state_t const & s, beam_state_t const & t) { return s.turn < t.turn; }  // weak
bool operator > (beam_state_t const & s, beam_state_t const & t) { return s.turn > t.turn; }

/// ビームサーチでUFOの軌道を構成します。
///
/// @param[in] turn_limit 打ち切りターン数
beam_state_t beamsearch(Stage const & stage, bitset<Parameter::MaxHouseCount> const & delivered, UFO const & ufo, int turn_limit) {
    vector<beam_state_t> beam; {
        beam_state_t initial = {};
        initial.delivered = delivered;
        initial.turn = 0;
        initial.house = -1;  // the office
        beam.push_back(initial);
    }
    int house_count = stage.houses().count();
    constexpr int beam_width = 30;
    set<pair<uint64_t, int> > used;
    while (true) {
        used.clear();
        vector<beam_state_t> prev_beam;
        prev_beam.swap(beam);
        for (auto const & s : prev_beam) {
            auto key = make_pair(hash<bitset<Parameter::MaxHouseCount> >()(s.delivered), s.house);
            if (used.count(key)) continue;
            used.insert(key);
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
        if (beam.size() < beam_width) {
            sort(whole(beam));
        } else {
            partial_sort(beam.begin(), beam.begin() + beam_width, beam.end());
            beam.resize(beam_width);
        }
        if (int(beam.front().delivered.count()) == house_count) break;
        if (beam.front().turn > turn_limit) break;
    }
    return beam.front();
}

/// ビームサーチで事前構築された移動計画に沿って、いい感じにActionsを構成します。
void move_items_with_large_ufos_plan(Stage const & stage, Actions & actions, TargetManager & target, vector<vector<int> > const & large_path, bitset<Parameter::MaxHouseCount> const & delivered_by_large) {
    int house_count = stage.houses().count();
    array<int, Parameter::UFOCount> item_count;
    repeat (ufo_index, Parameter::UFOCount) {
        auto const & ufo = stage.ufos()[ufo_index];
        item_count[ufo_index] = ufo.itemCount();

        // 農場に接触しているなら自明に補給すべき
        if (item_count[ufo_index] < ufo.capacity() and Util::IsIntersect(ufo, stage.office())) {
            actions.add(Action::PickUp(ufo_index));
            item_count[ufo_index] = ufo.capacity();
        }
        // 大きいUFOに接触しているならたいてい補給した方がよい
        if (item_count[ufo_index] < ufo.capacity() and ufo.type() == UFOType_Small) {
            repeat (large_ufo_index, Parameter::LargeUFOCount) {
                auto const & large_ufo = stage.ufos()[large_ufo_index];
                if (Util::IsIntersect(ufo, large_ufo)) {
                    actions.add(Action::Pass(large_ufo_index, ufo_index));
                    int delta = min(item_count[large_ufo_index], ufo.capacity() - item_count[ufo_index]);
                    item_count[ufo_index] += delta;
                    item_count[large_ufo_index] -= delta;
                }
            }
        }

        // アイテムないならロックを手放す
        if (item_count[ufo_index] == 0) {
            if (target.is_targetting(ufo_index)) {
                target.unlink_ufo(ufo_index);
            }

        } else {
            // 目標の家に着いたら配達
            if (target.is_targetting(ufo_index)) {
                int house_index = target.from_ufo(ufo_index);
                auto const & house = stage.houses()[house_index];
                if (Util::IsIntersect(ufo, house)) {
                    item_count[ufo_index] -= 1;
                    actions.add(Action::Deliver(ufo_index, house_index));
                    target.deliver_house(house_index);
                }
            }

            // 目標がないなら、大きなUFOは(存在するなら)計画に沿って移動
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

            // 目標がないなら、とりあえず近い家を宣言する
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

            // 目標がないなら、他のUFOから貪欲に奪う
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

        // アイテムがないなら最寄りの補給可能場所へ
        if (ufo.itemCount() == 0) {
            Vector2 pos = stage.office().pos();
            if (ufo.type() == UFOType_Small) {
                repeat (large_ufo_index, Parameter::LargeUFOCount) {
                    auto const & large_ufo = stage.ufos()[large_ufo_index];
                    if (large_ufo.itemCount() >= 5) {
                        if (ufo.pos().dist(large_ufo.pos()) < ufo.pos().dist(pos)) {
                            pos = large_ufo.pos();  // NOTE: The prediction is required.
                        }
                    }
                }
            }
            target_positions.add(pos);

        // 家へ向かう
        } else {
            int house_index = target.from_ufo(ufo_index);
            if (house_index != TargetManager::NONE) {
                target_positions.add(stage.houses()[house_index].pos());

            // 暇なら待機
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
#ifdef LOCAL
int current_stage = -1;
#endif

}

//------------------------------------------------------------------------------
/// 各ステージ開始時に呼び出されます。
///
/// ここで、各ステージに対して初期処理を行うことができます。
///
/// @param[in] stage 現在のステージ。
void Answer::init(Stage const & a_stage) {
    using namespace Solver;

    result.clear();
#ifdef LOCAL
    current_stage += 1;
#endif
detect_towns(a_stage.houses());

    for (int turn_limit : { 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140 }) {
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

#ifdef DEBUG
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

#ifdef DEBUG
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

#ifdef LOCAL
    const char *green = "\x1b[32m";
    const char *bright_green = "\x1b[32;1m";
    const char *bright_yellow = "\x1b[33;1m";
    const char *yellow = "\x1b[33m";
    const char *bright_red = "\x1b[31;1m";
    const char *red = "\x1b[31m";
    const char *color =
        int(result.size()) <= 50 ? green :
        int(result.size()) <= 70 ? bright_green :
        int(result.size()) < 100 ? bright_yellow :
        int(result.size()) < 130 ? yellow :
        int(result.size()) < 150 ? bright_red :
        red;
    fprintf(stderr, "%5d |%s%5d\x1b[0m\n", current_stage, color, int(result.size()));
#endif
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
    using namespace Solver;
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
    using namespace Solver;
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
