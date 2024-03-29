﻿//------------------------------------------------------------------------------
/// @file
/// @author   ハル研究所プログラミングコンテスト実行委員会
///
/// @copyright  Copyright (c) 2017 HAL Laboratory, Inc.
/// @attention  このファイルの利用は、同梱のREADMEにある
///             利用条件に従ってください。
//------------------------------------------------------------------------------

#pragma once

#include "Answer.hpp"
#include "Recorder.hpp"
#include "Random.hpp"
#include "Timer.hpp"

namespace hpc {

//------------------------------------------------------------------------------
/// ゲーム全体。
class Game
{
public:
    Game(RandomSeed aSeed);
    void changeSeed(RandomSeed aSeed);         ///< シード値を変更します。
    void run(Answer& aAnswer);                 ///< ゲームを実行します。
    const Recorder& recorder()const;           ///< ログ記録器を取得します。
    const Timer& timer()const;                 ///< タイマーを取得します。
private:
    Random mRandom;                            ///< 乱数生成器
    Recorder mRecorder;                        ///< ログ記録器
    Timer mTimer;                              ///< タイマー
};

} // namespace
// EOF
