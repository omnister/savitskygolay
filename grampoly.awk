
awk '
  BEGIN {
    print -2, -3/(-6+24+17)
    print -1, 12/(-6+24+17)
    print  0, 17/(-6+24+17)
    print 1, 12/(-6+24+17)
    print 2, -3/(-6+24+17)
    print -3, -2/(-4+6+12+7) 
    print -2, 3/(-4+6+12+7) 
    print -1, 6/(-4+6+12+7) 
    print 0, 7/(-4+6+12+7) 
    print 1, 6/(-4+6+12+7) 
    print 2, 3/(-4+6+12+7) 
    print 3, -2/(-4+6+12+7)2
    print -4, -21/(-42+28+78+108+59)
    print -3, 14/(-42+28+78+108+59)
    print -2, 39/(-42+28+78+108+59)
    print -1, 54/(-42+28+78+108+59)
    print 0, 59/(-42+28+78+108+59)
    print 1, 54/( -42+28+78+108+59)
    print 2, 39/(-42+28+78+108+59)
    print 3, 14/(-42+28+78+108+59)
    print 4, -21/(-42+28+78+108+59)
    print -5, -36/(-72+18+88+138+168+89)
    print -4, 9/(-72+18+88+138+168+89)
    print -3, 44/(-72+18+88+138+168+89)
    print -2, 69/(-72+18+88+138+168+89)
    print -1, 84/(-72+18+88+138+168+89)
    print -0, 89/(-72+18+88+138+168+89)
    print 1, 84/(-72+18+88+138+168+89)
    print 2, 69/(-72+18+88+138+168+89)
    print 3, 44/(-72+18+88+138+168+89)
    print 4, 9/(-72+18+88+138+168+89)
    print 5, -36/(-72+18+88+138+168+89)
  }
'

