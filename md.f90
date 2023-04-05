      PROGRAM md
      IMPLICIT NONE
      
!     実数型の変数：E(全エネルギー)，area(サイズ), ke(運動エネ), kecum(運動エネの和), pe(位置エネ), pecum(位置エネの和), t(時刻), vcum(速度の和), virial(圧力), Lx(縦), Ly(横),
!            ax(加速度x成分), ay(加速度y成分), dt(時刻の微小変化量), dt2(dt^2), vx(速度x成分), vy(速度y成分), x(x座標), y(y座標),?_old(1ステップ前の？)
      REAL*8 E, area, ke, kecum, pe, pecum, t, vcum, virial, Lx, Ly,&
            ax(36), ay(36),ax_old(36),ay_old(36), dt, dt2, vx(36), vy(36), x(36), y(36), dx_old(36), dy_old(36),&
             fxij, fyij

!     整数型の変数：IC1_,IC2_(描画ウィンドウ)，ncum(データ積算回数)
      INTEGER IC1_, IC2_, N, ncum
      
!     文字型の変数：flag(描画処理)
      CHARACTER*80 flag
      
!     論理型の変数：Ltmp1(繰り返し処理)
      LOGICAL Ltmp1_
      
      CALL initial()           !初期化する
      CALL set_up_windows()    !ウィンドウの設定
      CALL accel()             !加速度の計算
      E = ke + pe              !エネルギー＝運動エネ＋位置エネ
      flag = ''
      Ltmp1_ = .TRUE.
      DO WHILE(Ltmp1_)
         CALL show_positions()  !粒子の軌跡を描画
         CALL Verlet()          !速度型ベルレ法による速度計算
         CALL show_output()     !エネルギーの描画
         Ltmp1_ = .NOT.(flag .EQ. 'stop')
      END DO
      
      !ウィンドウを閉じる
      CALL gclose(IC1_)
      CALL gclose(IC2_)
      
      CONTAINS
      
      !初期条件の設定を行う
      SUBROUTINE initial()
         IMPLICIT NONE
         REAL*8 DATA_(64)
            
         INTEGER IO1_, i
         CHARACTER*80 file, heading, response, start
         LOGICAL Ltmp1_
         INTEGER DP_
         DATA DATA_/&
            1.09D0,0.98D0,-0.33D0,0.78D0,3.12D0,5.25D0,0.12D0,-1.19D0,&
            0.08D0,2.38D0,-0.08D0,-0.10D0,0.54D0,4.08D0,-1.94D0,-0.56D0,&
            2.52D0,4.39D0,0.75D0,0.34D0,3.03D0,2.94D0,1.70D0,-1.08D0,&
            4.25D0,3.01D0,0.84D0,0.47D0,0.89D0,3.11D0,-1.04D0,0.06D0,&
            2.76D0,0.31D0,1.64D0,1.36D0,3.14D0,1.91D0,0.38D0,-1.24D0,&
            0.23D0,5.71D0,-1.58D0,0.55D0,1.91D0,2.46D0,-1.55D0,-0.16D0,&
            4.77D0,0.96D0,-0.23D0,-0.83D0,5.10D0,4.63D0,-0.31D0,0.65D0,&
            4.97D0,5.88D0,1.18D0,1.48D0,3.90D0,0.20D0,0.46D0,-0.51D0/
         DATA DP_/1/, IO1_/11/
         dt = 0.01D0
         dt2 = dt**2
         response = ''
         Ltmp1_ = .TRUE.
         DO WHILE(Ltmp1_)
         !初期座標，初期速度について上記のデータか別のデータファイルから読み取るか選択
            WRITE(*,'(A,$)') 'read data statements (d) or file (f)? '
            READ(*,'(A)') start
         
            !上記のデータを初期データとするとき
            IF((start .EQ. 'd') .OR. (start .EQ. 'D')) THEN
               response = 'ok'
               
               !Lx×LyのサイズでN個の粒子を動かす
               !N個の粒子の初期座標，初期速度を読み取る
               N = 16
               Lx = 6
               Ly = 6
               DO i = 1, N
                  x(i) = DATA_(DP_)
                  DP_ = DP_ + 1
                  y(i) = DATA_(DP_)
                  DP_ = DP_ + 1
                  vx(i) = DATA_(DP_)
                  DP_ = DP_ + 1
                  vy(i) = DATA_(DP_)
                  DP_ = DP_ + 1
               END DO
            
            !別のデータファイルから読み込むとき
            !読み取りファイルは粒子数→xサイズ→yサイズ→座標：x，y→速度：vx，vyの順であること
            ELSE IF((start .EQ. 'f') .OR. (start .EQ. 'F')) THEN
               response = 'ok'
               WRITE(*,'(A,$)') 'file name = '
               READ(*,'(A)') file
               OPEN(IO1_, file = file)
               READ(IO1_,*) N
               READ(IO1_,*) Lx
               READ(IO1_,*) Ly
               READ(IO1_,'(A)') heading
               DO i = 1, N
                  READ(IO1_,*) x(i),y(i)
               END DO
               READ(IO1_,'(A)') heading
               DO i = 1, N
                  READ(IO1_,*) vx(i),vy(i)
               END DO
               CLOSE(IO1_)
            
            !d(D),f(F)以外の文字が入力されたら再度入力を求める
            ELSE
               WRITE(*,*)
               WRITE(*,*) 'd or f are the only acceptable responses.'
            END IF
            Ltmp1_ = response .NE. 'ok'
         END DO
         
   !     ke: 運動エネルギー
         ke = 0.0D0
         
   !     全ての粒子の1ステップ前のdx, dy, 運動エネルギーを更新
         DO i = 1, N
            dx_old(i) = vx(i)*dt
            dy_old(i) = vy(i)*dt
            ke = ke + vx(i)**2 + vy(i)**2
         END DO
         ke = 0.5D0*ke
         area = Lx*Ly
   !     t: 時刻
         t = 0.0D0
   !     和の初期化
         kecum = 0.0D0
         pecum = 0.0D0
         vcum = 0.0D0
      END SUBROUTINE initial
      
!     ウィンドウの設定
      SUBROUTINE set_up_windows()
         IMPLICIT NONE
   !     数値の出力ウィンドウ
         CALL gopen(600, 40, IC1_)
         CALL gsetbgcolor(IC1_,'Blue'//CHAR(0))
         CALL gclr(IC1_)
         CALL newwindow(IC1_, 0.0, 4.0, 60.0, 0.0)
         CALL headings()
   !     粒子の軌跡ウィンドウ
         CALL gopen(600, 600, IC2_)
         CALL newwindow(IC2_, -0.1*Real(Lx),-0.1*Real(Ly), &
                        Real(Lx)*1.1, Real(Ly)*1.1)
         CALL newpencolor(IC2_, 0)
         CALL gsetbgcolor(IC2_,'white'//CHAR(0))
         CALL gclr(IC2_)
         CALL drawrect(IC2_, 0.0, 0.0, Real(Lx), REAL(Ly))
      END SUBROUTINE set_up_windows
      
!     見出し
      SUBROUTINE headings()
         IMPLICIT NONE
         CALL drawstr(IC1_, 1.0, 1.5, 16.0,'time steps ',0.0,11)
         CALL drawstr(IC1_, 12.0, 1.5, 16.0,'time ',0.0,5)
         CALL drawstr(IC1_, 24.0, 1.5, 16.0,'energy ',0.0,7)
         CALL drawstr(IC1_, 36.0, 1.5, 16.0,'<T> ',0.0,4)
         CALL drawstr(IC1_, 48.0, 1.5, 16.0,'<P> ',0.0,4)
      END SUBROUTINE headings
      
!     速度型ベルレ法
      SUBROUTINE Verlet()
         IMPLICIT NONE
         REAL*8 xnew, ynew, dxnew, dynew
         INTEGER i
   !
         ke = 0.0d0
   !
         DO i = 1, N
   !     最新の座標と現在の座標の差を最初に求める
            dxnew =  dx_old(i) + ax(i)*dt2
            dynew =  dy_old(i) + ay(i)*dt2
   !     最新の座標を求める
            xnew = x(i) + dxnew
            ynew = y(i) + dynew
   !     周期境界条件を適用し，新しい座標を確定する
            x(i) = pbc(xnew, Lx)
            y(i) = pbc(ynew, Ly)
   !     最新の速度を現在の速度と加速度から求める
            vx(i) = vx(i) + (ax(i) + ax_old(i))*0.5d0*dt
            vy(i) = vy(i) + (ay(i) + ay_old(i))*0.5d0*dt
   !     新しい運動エネルギー
            ke = ke + vx(i)**2 + vy(i)**2
   !     最新の座標と現在の座標の差を保存する
            dx_old(i) = dxnew
            dy_old(i) = dynew
         END DO
   !     新しい加速度を求める
         ax_old = ax
         ay_old = ay
         CALL accel()
   !
         ke = 0.5D0*ke
         t = t + dt
      END SUBROUTINE Verlet
      
!     周期境界条件:容器から出た粒子は同じ速度で早退する位置から出てくる
      FUNCTION pbc(pos,L)
         IMPLICIT NONE
         REAL*8 pbc, pos, L
         IF(pos .LT. 0) THEN
            pbc = pos + L
         ELSE IF(pos .GT. L) THEN
            pbc = pos - L
         ELSE
            pbc = pos
         END IF
      END FUNCTION pbc
      
!     加速度を計算
      SUBROUTINE accel()
         IMPLICIT NONE
         REAL*8 dx, dy, pot
         INTEGER i, j
         DO i = 1, N
            ax(i) = 0.0d0
            ay(i) = 0.0d0
         END DO
         pe = 0.0d0
         virial = 0.0d0
         DO i = 1, N - 1
            DO j = i + 1, N
               dx = separation(x(i) - x(j),Lx)
               dy = separation(y(i) - y(j),Ly)
   !           換算単位系では 質量 = 1 なので 加速度 = 力 である
               CALL force(dx,dy,fxij,fyij,pot)
               ax(i) = ax(i) + fxij
               ay(i) = ay(i) + fyij
               ax(j) = ax(j) - fxij
               ay(j) = ay(j) - fyij
               pe = pe + pot
               virial = virial + dx*fxij + dy*fyij
            END DO
         END DO
      END SUBROUTINE accel
         
   !     力を計算
      SUBROUTINE force(dx,dy,fx,fy,pot)
         IMPLICIT NONE
         REAL*8 f_over_r, r2, rm2, rm6, dx, dy, fx, fy, pot
         r2 = dx*dx + dy*dy
         rm2 = 1/r2
         rm6 = rm2*rm2*rm2
         f_over_r = 24*rm6*(2*rm6 - 1)*rm2
         fx = f_over_r*dx
         fy = f_over_r*dy
         pot = 4.0d0*(rm6*rm6 - rm6)
      END SUBROUTINE force
         
   !     
      FUNCTION separation(ds,L)
         IMPLICIT NONE
         REAL*8 separation, ds, L
         IF(ds .GT. 0.5D0*L) THEN
            separation = ds - L
         ELSE IF(ds .LT. -0.5D0*L) THEN
            separation = ds + L
         ELSE
            separation = ds
         END IF
      END FUNCTION separation
      
!     粒子の軌跡の描画：入力されたキーによって描画方法を変える
      SUBROUTINE show_positions()
         IMPLICIT NONE
         INTEGER i, k
         CALL ggetch(k)
   !     再描画して軌跡位置を描画（軌跡を消去）　'r': 114
         flag = ''
         IF(k .EQ. 114) THEN
            CALL gclr(IC2_)
            CALL drawrect(IC2_, 0.0, 0.0, Real(Lx), REAL(Ly))
   !     終了　's': 115
         ELSE IF(k .EQ. 115) THEN
            flag = 'stop'
   !     粒子位置を描画しない　'n': 110
         ELSE IF(k .EQ. 110) THEN
            flag = 'no_show'
         END IF
         CALL gclr(IC1_)
         CALL headings()
         IF(flag .NE. 'no_show') THEN
            DO i = 1, N
               CALL fillcirc(IC2_, REAL(x(i)), REAL(y(i)), 0.07, 0.07  )
            END DO
         END IF
      END SUBROUTINE show_positions
      
!     エネルギーなどの数値を描画
      SUBROUTINE show_output()
         IMPLICIT NONE
         REAL*8 E, mean_ke, p
         ncum = ncum + 1
   !     データ積算回数
         CALL drawnum(IC1_, 1.0, 3.5, 16.0, REAL(ncum),0.0, -1)
   !     時間
         CALL drawnum(IC1_, 12.0, 3.5, 16.0, REAL(t),0.0, 2)
   !     全エネルギー = 運動エネ＋位置エネ
         E = ke + pe
         CALL drawnum(IC1_, 24.0, 3.5, 16.0, REAL(E),0.0, 4)
         kecum = kecum + ke
         vcum = vcum + virial
         mean_ke = kecum/ncum
   !     さらに N で割る必要がある
         p = mean_ke + (0.5D0*vcum)/ncum
   !     平均圧力
         p = p/area
   !     平均の運動学的温度
         CALL drawnum(IC1_, 36.0, 3.5, 16.0, Real(mean_ke/N),0.0, 4)
   !     平均圧力
         CALL drawnum(IC1_, 48.0, 3.5, 16.0, Real(p),0.0, 4)
      END SUBROUTINE show_output
      END