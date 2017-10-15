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
#include <random>
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

constexpr float eps = 1e-4;

minstd_rand gen;

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
bool operator < (town_t const & a, town_t const & b) {
    if (a.center == b.center) {
        return a.house_indices < b.house_indices;
    }
    return make_pair(a.center.y, a.center.x) < make_pair(b.center.y, b.center.x);
}

vector<town_t> reconstruct_towns_from_centers(vector<Vector2> const & town_centers, int radius, Houses const & houses) {
    vector<town_t> towns;
    for (auto town_center : town_centers) {
        town_t town = {};
        town.center = town_center;
        repeat (house_index, houses.count()) {
            auto const & house = houses[house_index];
            if (house.pos().dist(town.center) <= radius) {
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
    return towns;
}

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
    vector<town_t> towns = reconstruct_towns_from_centers(town_centers, StageParameter::TownRadius + 3, houses); // 3 は余裕

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

vector<Vector2> get_town_centers(vector<town_t> const & towns) {
    vector<Vector2> poss;
    for (auto const & town : towns) {
        poss.push_back(town.center);
    }
    return poss;
}

/// 2点間の移動に要するターン数を計算します。
///
/// @note 改善余地あり。半径とかをまったく考慮していない。
int turns_to_move(Vector2 const & from, Vector2 const & to, int max_speed) {
    return ceil(from.dist(to) / max_speed);
}
pair<int, Vector2> simulate_move(Vector2 const & from, Vector2 const & to, int sum_radius, int max_speed) {
    Vector2 pos = from;
    int turn = 0;
    for (; ; ++ turn) {
        Vector2 dir = to - pos;
        float length = dir.length();
        if (length < sum_radius - eps) {  // 保守的のためeps
            break;
        }
        if (length > max_speed) {
            dir.unitAssign(max_speed - eps);
        }
        pos = pos + dir;
    }
    return { turn, pos };
}

/// 配達計画から実際の経路や補給を構成する
pair<vector<Vector2>, vector<pair<int, Action> > > get_delivering_path(vector<int> const & delivering_plan, Stage const & stage, int ufo_index, vector<vector<Vector2> > const & large_ufo_path) {
    auto const & ufo = stage.ufos()[ufo_index];
    vector<Vector2> positions;
    vector<pair<int, Action> > actions;
    Vector2 pos = stage.office().pos();
    auto move_to = [&](Vector2 const & target_pos, int target_radius) {
        while (true) {
            Vector2 dir = target_pos - pos;
            float length = dir.length();
            if (length < ufo.radius() + target_radius - eps) break;
            if (length > ufo.maxSpeed()) dir.unitAssign(ufo.maxSpeed() - eps);
            pos = pos + dir;
            positions.push_back(pos);
        }
    };
    // 大きいUFOの位置を取得
    auto get_large_ufo_pos = [&](int large_ufo_index, int turn) {
        // NOTE: (他を実装するとき)ここで配列外だった場合に注意する
        if (large_ufo_path.size() <= large_ufo_index) return stage.office().pos();
        auto const & large_path = large_ufo_path[large_ufo_index];
        return
            large_path.empty() or turn <= 0 ? stage.office().pos() :
            large_path.size() <= turn - 1 ? large_path.back() :
            large_path[turn - 1];
    };
    auto chase_large_ufo = [&](int large_ufo_index) {
        while (true) {
            Vector2 dir = get_large_ufo_pos(large_ufo_index, positions.size()) - pos;
            float length = dir.length();
            if (length < ufo.radius() + Parameter::LargeUFORadius - eps) break;
            if (length > ufo.maxSpeed()) dir.unitAssign(ufo.maxSpeed() - eps);
            pos = pos + dir;
            positions.push_back(pos);
        }
    };

    // 初手で積み込みは確定
    actions.emplace_back(0, Action::PickUp(ufo_index));
    int item_count = ufo.capacity();
    // 前から見ていく
    for (int house_index : delivering_plan) {
        auto const & house = stage.houses()[house_index];

        // 補充が必要な場合
        if (item_count == 0) {
            int target_index = -1;
            if (ufo.type() == UFOType_Small) {
                int required_turn = turns_to_move(pos, stage.office().pos(), Parameter::SmallUFOMaxSpeed) + turns_to_move(stage.office().pos(), house.pos(), Parameter::SmallUFOMaxSpeed);
                repeat (large_ufo_index, Parameter::LargeUFOCount) {
                    // 移動先の予測
                    Vector2 large_pos = get_large_ufo_pos(large_ufo_index, positions.size());
                    int large_turn = turns_to_move(pos, large_pos, Parameter::SmallUFOMaxSpeed) + turns_to_move(large_pos, house.pos(), Parameter::SmallUFOMaxSpeed);
                    if (large_turn < required_turn) {
                        target_index = large_ufo_index;
                    }
                }
            }
            if (target_index == -1) {
                move_to(stage.office().pos(), Parameter::OfficeRadius);
                actions.emplace_back(positions.size(), Action::PickUp(ufo_index));
            } else {
                chase_large_ufo(target_index);
                actions.emplace_back(positions.size(), Action::Pass(target_index, ufo_index));
            }
            item_count = ufo.capacity();
        }

        // 配りに行く
        move_to(house.pos(), Parameter::HouseRadius);
        actions.emplace_back(positions.size(), Action::Deliver(ufo_index, house_index));
        item_count -= 1;
    }
    if (not actions.empty()) {
        while (positions.size() <= actions.back().first) {
            positions.push_back(pos);
        }
    }
    assert (not positions.empty());
    return { positions, actions };
}

/// 街に所属しない家の一覧を取得
vector<int> get_countryside_house_indices(int house_count, vector<town_t> const & towns) {
    vector<int> xs(house_count);
    iota(whole(xs), 0);
    for (auto const & town : towns) {
        for (int house_index : town.house_indices) {
            xs[house_index] = -1;
        }
    }
    xs.erase(remove(whole(xs), -1), xs.end());
    return xs;
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
    vector<town_t> towns = detect_towns(a_stage.houses());
#ifdef LOCAL
    current_stage += 1;
#endif

    towns = reconstruct_towns_from_centers(get_town_centers(towns), StageParameter::TownRadius * 2, a_stage.houses());
    sort(whole(towns));

    vector<vector<int> > delivering_plan(Parameter::UFOCount);
    repeat (house_index, a_stage.houses().count()) {
        delivering_plan[house_index % Parameter::UFOCount].push_back(house_index);
    }
    vector<vector<Vector2> > positions_list;
    vector<vector<pair<int, Action> > > actions_list;
    repeat (ufo_index, Parameter::UFOCount) {
        vector<Vector2> positions;
        vector<pair<int, Action> > actions;
        tie(positions, actions) = get_delivering_path(delivering_plan[ufo_index], a_stage, ufo_index, positions_list);
        positions_list.push_back(positions);
        actions_list.push_back(actions);
    }
    int result_turn = 0;
    repeat (ufo_index, Parameter::UFOCount) {
        setmax(result_turn, int(positions_list[ufo_index].size()));
    }
    result.resize(result_turn);
    repeat (ufo_index, Parameter::UFOCount) {
        for (auto it : actions_list[ufo_index]) {
            result[it.first].actions.add(it.second);
        }
        repeat (turn, result_turn) {
            result[turn].target_positions.add(
                turn < positions_list[ufo_index].size()
                ? positions_list[ufo_index][turn]
                : result[turn - 1].target_positions[ufo_index]);
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
// if (stage.turn() >= result.size()) return;
    assert (stage.turn() < result.size());
    actions = result[stage.turn()].actions;
    for (auto action : actions) {
        switch (action.type()) {
            case ActionType_PickUp:
                assert (Util::IsIntersect(stage.ufos()[action.ufoIndex()], stage.office()));
                break;
            case ActionType_Pass:
                assert (stage.ufos()[action.srcUFOIndex()].itemCount() >= stage.ufos()[action.dstUFOIndex()].capacity());
                assert (Util::IsIntersect(stage.ufos()[action.srcUFOIndex()], stage.ufos()[action.dstUFOIndex()]));
                break;
            case ActionType_Deliver:
                assert (Util::IsIntersect(stage.ufos()[action.ufoIndex()], stage.houses()[action.houseIndex()]));
                break;
            default:
                break;
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
    using namespace Solver;
// if (stage.turn() >= result.size()) { target_positions = result.back().target_positions; return; }
    target_positions = result[stage.turn()].target_positions;
    if (stage.turn() >= 1) {
        repeat (ufo_index, Parameter::UFOCount) {
            assert (stage.ufos()[ufo_index].pos().dist(result[stage.turn() - 1].target_positions[ufo_index]) < eps);
        }
    }
#ifdef DEBUG
repeat (house_index, stage.houses().count()) cerr << stage.houses()[house_index].delivered();
cerr << endl;
#endif
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
