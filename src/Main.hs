{-# LANGUAGE ScopedTypeVariables, RankNTypes , FlexibleInstances #-}
module Main where

import qualified Data.Map.Strict as Map

import SolutionDomain
import CalculusUtils
import Data.List
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as M
import Data.Maybe

--{-
import Control.Monad.Par as Par
import Control.Monad
import Control.Monad.Par.IO as ParIO
import Control.Monad.Trans (liftIO)
import qualified Graphics.Gnuplot.Graph as Graph
import qualified Graphics.Gnuplot.Frame as Frame
import qualified Graphics.Gnuplot.Frame.OptionSet as Opts
import qualified Graphics.Gnuplot.Frame.Option as Opt
import qualified Graphics.Gnuplot.Advanced as GP
import qualified Graphics.Gnuplot.Plot.TwoDimensional as Plot2D
import qualified Graphics.Gnuplot.Graph.TwoDimensional as Graph2D
import Graphics.Gnuplot.Terminal.PNG as PNG
import Control.DeepSeq
import GHC.Generics (Generic)
import qualified Control.DeepSeq

instance Par.NFData Property
instance Par.NFData AdjNode



defltOpts :: Graph.C graph => Opts.T graph
defltOpts = Opts.key False Opts.deflt

plotDomain grid scale=
    let xsize = length (head grid)-1
        ysize = length  grid  -1
        printSize = show (ysize * 5) ++ "," ++ show (xsize * 5)
    in Frame.cons
            (   Opts.add ( Opt.custom "terminal pngcairo size" printSize) [printSize] $
                scale $
                Opts.sizeRatio (fromIntegral xsize / fromIntegral ysize)  defltOpts
            ) $
        Plot2D.function Graph2D.image
        (liftM2 (,)  [0..ysize] [0..xsize])
            $ \(x,y) -> (grid!!x)!!y

plotDomainLog grid = plotDomain grid $ Opts.add ( Opt.custom "logscale cb"  "") []

plotDomainLinear grid = plotDomain grid id
---}

heatConductivityK:: Double
heatConductivityK = 0.025

orthogonalSides:: Side ->[Side]
orthogonalSides side = case side of
    East -> [North, South, Top, Bottom]
    West -> [North, South, Top, Bottom]
    North -> [East, West, Top,Bottom]
    South -> [East, West, Top,Bottom]
    Top -> [East,West,North,South]
    Bottom -> [East,West,North,South]
    Now -> [Now,Prev]
    Prev -> [Now,Prev]
    Center -> [Center]

addTerm:: [Term a]->Term a->[Term a]
addTerm terms term = terms ++ [term]

addTerms:: [Term a]-> [Term a]-> [Term a]
addTerms terms1 terms2 = case terms2 of
    ([x]) -> addTerm terms1 x
    (x:xs) -> addTerms (addTerm terms1 x) xs
    _ -> terms1

testEquation:: (Num a, Fractional a ) => Equation (Term a )
testEquation =
    Equation
        --[Unknown 2, Constant 4, Constant (-0.232)]
        [Unknown 2, SubExpression (Expression [Constant 2, Constant 2, Unknown 0.025]), Constant (-0.232)]
        ( addTerms [Unknown (-0.025),Unknown (-0.05)] [Constant 2, Constant 3 ] )
        U

--continuity::(Num a, Fractional a, RealFloat a)=> Reader (ValSet a) (Equation (Position-> [Term a]))
--continuity::Reader AdjGraph (Equation (AdjPos-> [Term Double]))
continuity env = let integrate = integUnknown env Density
        in Equation
            ((integrate Temporal [drho_dt env] >>= integrate Spatial)
              : integrateTerms integrate env ( divergenceWithProps [Density]) )
            [ const [Constant 0]]
            Density

--uMomentum::(Num a, Fractional a, RealFloat a)=> Reader (ValSet a) (Equation (Position-> [Term a]))
--uMomentum::Reader (ValSet Double) (Equation (Position-> [Term Double]))
uMomentum env = let integrate = integUnknown env U
        in Equation
            (  [ integrate Temporal [drhou_dt env] >>= integrate Spatial ]
                ++ integrateTerms integrate env (divergenceWithProps [Density, U])
                 ++ [integrate Spatial [dp_dx env] >>= integrate Temporal ]
            )
            ( concatMap (integrateTerms integrate env)
                [  divGrad [Mew,U] 1
                  ,divergence [ dmewu_dx, dmewv_dx, dmeww_dx ]
                  ,map (\d-> ddf_ d (-2/3) X) (divergenceWithProps [Mew])
                ]
            )
            U

--vMomentum::Reader (ValSet Double) (Equation (Position-> [Term Double]))
vMomentum env = let integrate = integUnknown env V
        in Equation
            ([ integrate Temporal [ drhov_dt env] >>= integrate Spatial ]
                ++ integrateTerms integrate env (divergenceWithProps [Density, V])
                ++ [integrate Spatial [ dp_dy env] >>= integrate Temporal ]
            )
            ( concatMap (integrateTerms integrate env)
                [ divGrad [Mew,V] 1
                  ,divergence [ dmewu_dy, dmewv_dy, dmeww_dy ]
                  ,map (\d-> ddf_ d (-2/3) Y) (divergenceWithProps [Mew])
                ]
            )
            V

--wMomentum::Reader (ValSet Double) (Equation (Position-> [Term Double]))
wMomentum env =  let integrate = integUnknown env W
        in Equation
            ([ integrate Temporal [ drhow_dt env] >>= integrate Spatial ]
                ++ integrateTerms integrate env (divergenceWithProps [Density , W] )
                ++[integrate Spatial [ dp_dz env] >>= integrate Temporal ]
            )
            ( concatMap (integrateTerms integrate env)
                [ divGrad [Mew,W] 1
                  ,divergence [ dmewu_dz, dmewv_dz, dmeww_dz ]
                  ,map (\d-> ddf_ d (-2/3) Z) (divergenceWithProps [Mew])
                ]
            )
            W

dye env = let integrate = integUnknown env Dye
        in Equation
            ([ integrate Temporal [ drhodye_dt env] >>= integrate Spatial ]
                ++ integrateTerms integrate env (divergenceWithProps [Density , Dye] )
            )
            ( concatMap (integrateTerms integrate env)
                [ divGrad [Mew,Dye] 0.001  ]
            )
            Dye

--energy::Reader (ValSet Double) (Equation (Position-> [Term Double]))
energy env = let integrate = integUnknown env Temperature
        in Equation
            ([ integrate Temporal [ drhoT_dt env] >>= integrate Spatial ]
                ++ integrateTerms integrate env (divergenceWithPropsFactor [Density,Temperature] specificHeatCv )
                ++ integrateTerms integrate env (divergenceWithProps [Pressure] )
            )
            ( concatMap (integrateTerms integrate env)
                [ divGrad [Temperature] heatConductivityK
                  , squareDerivative [Mew,U] (8/3) X
                  , squareDerivative [Mew,V] (8/3) Y
                  , squareDerivative [Mew,W] (8/3) Z
                  , squareDerivative [Mew,U] 1 Y
                  , squareDerivative [Mew,V] 1 X
                  , squareDerivative [Mew,U] 1 Z
                  , squareDerivative [Mew,W] 1 X
                  , squareDerivative [Mew,V] 1 Z
                  , squareDerivative [Mew,W] 1 Y
                  , pairedMultipliedDerivatives [Mew,U][V] X Y
                  , pairedMultipliedDerivatives [Mew,U][W] X Z
                  , pairedMultipliedDerivatives [Mew,W][V] Z Y
                ]
            )
            Temperature

getDiscEqInstance :: forall a t.Equation (t -> [a]) -> t -> Equation a
getDiscEqInstance (Equation l r up) pos = Equation (concatMap (\t -> t pos) $! l) (concatMap (\t -> t pos) $! r) up

calcSteps_adj = [
    (continuity,True) --run on all nodes
    ,(uMomentum,False)
    ,(vMomentum,False)
    ,(wMomentum,False)
    ,(energy,  True)  --run on all nodes
    ,(dye, False)
    ]

--setNodeProperty :: AdjNode -> Property -> Double -> IO AdjNode
setNodeProperty (AdjNode f e p n o a i) prp value tf = do
    --putStrLn $ "wrote value at node: " ++ show i ++ " "++ show o ++ " " ++ show value ++ " " ++ show prp
    --return $ AdjNode f e (p V.// [(propIndex "from setNodeProp" prp ,Just value)] ) n o a i
    --{-
    newp <- tf p
    M.write newp (propIndex "from setNodeProp" prp) $ Just value
    --putStrLn $ "wrote value at node: " ++ show i ++ " "++ show o ++ " " ++ show value ++ " " ++ show prp
    done <- V.unsafeFreeze newp
    return $ AdjNode f e done n o a i
    ---}

runSingleStep :: AdjGraph -> IO AdjGraph
runSingleStep adjgraph = foldM'
    (\prev next ->
        let stepTime = next == 0
        in do
            putStrLn $ "starting predictor step with stepTime: " ++ show stepTime
            predicted <- applyDiffEq_adj calcSteps_adj stepTime prev
            putStrLn $ "starting corrector step "
            correctPrevTime predicted calcSteps_adj
            --return predicted
    )
    adjgraph [0..2]

referencesSelf node self = 
    let n = neighbors node
        neighs = map sideIndex [East,West,North,South] 
    in foldl' (\prev next -> prev || self == ( fromJust $ n V.! next)) False neighs  

applyDiffEq_adj :: [(AdjGraph -> Equation ((Int, Int) -> [Term Double]), Bool)]
                         -> Bool -> AdjGraph -> IO AdjGraph
applyDiffEq_adj equations stepTime (AdjGraph aln cmt cat )=
    let new_cmt = if stepTime then pushUpTime cmt else cmt
        new_cat = if stepTime then cat + 1 else cat
        source_grid = grids aln V.! cmt
        env = AdjGraph aln cmt cat
    in do
        putStrLn $ "predictor found target grid with index: " ++ show new_cmt
        sg <- V.thaw $ grids aln V.! new_cmt
        _ <- ParIO.runParIO $ parMapM (\node -> liftIO $
                let j = index node
                    pos = (j,cmt) -- cmt is used here because regardless of if you save it here or next time slot the source data is always this slot
                in do
                    --putStrLn $ "predictor working on node with index: " ++ show j ++ " "++ (show $ origPos node)
                    foldM' (\prev (next_eq, alwaysRun) ->
                            let discEquation = force $ getDiscEqInstance (next_eq env) pos
                                solvedProperty = force $ unknownProperty discEquation
                                newValue = force $ solveUnknown discEquation pos env
                            in do
                                newNode <- setNodeProperty prev solvedProperty newValue V.thaw
                                M.write sg j (if (alwaysRun || active prev) -- && not (referencesSelf node j)  
                                    then newNode else prev)
                                --print "wrote results to node in predictor"
                                M.read sg j
                        ) node equations
            ) source_grid
        ng <- V.unsafeFreeze sg
        putStrLn $ "totalDensity in predictor: "
            ++ (show $ V.foldl (\prev (AdjNode _ _ p _ _ _ _) -> prev + (fromJust $ p V.! 3) ) 0.0 ng)
        return.force $ AdjGraph
            (Grids $ grids aln V.// [(new_cmt, ng)])
            new_cmt new_cat

correctPrevTime :: AdjGraph -> [(AdjGraph -> Equation ((Int, Int) -> [Term Double]), Bool)] -> IO AdjGraph
correctPrevTime (AdjGraph aln cmt cat ) equations =
    if cat > 1
        then
            let to_be_corrected_t = pushBackTime cmt
                concrete_t = pushBackTime to_be_corrected_t
                prev_grid = grids aln V.! to_be_corrected_t
                env = AdjGraph aln cmt cat
            in do
                putStrLn $ "corrector found target grid with index: " ++ show to_be_corrected_t
                sg <- V.unsafeThaw $ grids aln V.! to_be_corrected_t
                _ <- ParIO.runParIO $ parMapM (\node -> liftIO $
                        let j = index node
                        in do
                            --putStrLn $ "corrector working on node with index: " ++ show j
                            foldM' (\prev (next_eq, alwaysRun) ->
                                    let solvedProperty = unknownProperty (next_eq env)
                                        c = force $ prop_adj Nondirectional solvedProperty (j,cmt) Center env
                                        pp = force $ prop_adj Nondirectional solvedProperty (j,concrete_t) Center env
                                        newValue = force $ pp + (0.5 * ( c - pp))
                                    in do
                                        newNode <- setNodeProperty prev solvedProperty newValue V.unsafeThaw
                                        M.write sg j (if (alwaysRun || active prev) -- && not (referencesSelf node j)
                                            then newNode else prev)
                                        --print "wrote results to node in predictor"
                                        M.read sg j
                                ) node equations
                    ) prev_grid
                ng <- V.unsafeFreeze sg
                return.force $ AdjGraph
                    (Grids $ grids aln V.// [(to_be_corrected_t, ng)])
                    cmt cat
        else return.force $ AdjGraph aln cmt cat

supportCalcSteps = []

allPositionsCurrentTime env =
    map (\(idx,_)-> (idx, currModTime env) ) makeAllPositions_adj

prevGrid g = (grids $ allNodes g) V.! (pushBackTime $ currModTime g)

currGrid g = (grids $ allNodes g) V.! (currModTime g)

totalDensity env= foldl' (\p n -> p + prop_adj Nondirectional Density n Center env ) 0.0
    $ map (\x-> (x, pushBackTime $ currModTime env) ) [0.. (V.length $ prevGrid env )-1]

foldM' :: (Monad m) => (a -> b -> m a) -> a -> [b] -> m a
foldM' _ z [] = return z
foldM' f z (x:xs) = do
  z' <- f z x
  z' `seq` foldM' f z' xs

runTimeSteps_Print :: IO AdjGraph
runTimeSteps_Print =
    foldM'
        (\prev step -> do
            putStrLn $ "about to run step: " ++ show step
            next <-  runSingleStep prev
            {-
            -- use the most recent valset but one time step BACK because that will have been corrected. the most recent time level is NOT YET CORRECTED. it has JUST BEEN CALCULATED AND IS INTERMEDIATE
            mapM_
                (\nextProp ->
                    GP.plotAsync ( PNG.cons $ "c:\\temp\\"++ show nextProp ++ "\\"++show step ++".png")
                    $! plotDomainLinear $! adjGraphToGrid next  nextProp (1+maxPos Y) (quot (maxPos Z)  2)
                ) $ enumFrom Speed
            mapM_
                (\nextProp ->
                    GP.plotAsync ( PNG.cons $ "c:\\temp\\"++ show nextProp ++ " LOGSCALE" ++ "\\"++show step ++".png")
                    $! plotDomainLog $! adjGraphToGrid next  nextProp (1+maxPos Y) (quot (maxPos Z)  2)
                ) $ enumFrom Speed
            -}
            putStrLn $ show $ length (allPositionsCurrentTime next)
            putStrLn $ "step: " ++ show step
            putStrLn $ "timeLevel: " ++ show (currModTime next)
            putStrLn $ "timeLevelAbsolute: " ++ show (currAbsTime next)
            putStrLn $ "totalDensity: " ++ show (totalDensity next)
            putStrLn " "
            --putStrLn $ stringDomain U timeLevel (1 + maxPos Y) prev (quot (maxPos Z)  2)
            putStrLn "Forcing evaluation of graph..."
            return $ force next
        )
        initialGrid_adj
        [0..9999999]

{-
testTerms = [Unknown 2.4, Constant 1.2, Constant 3.112, Unknown (-0.21),  SubExpression (Expression [Constant 2, Constant 2, SubExpression (Expression [Unknown 0.33333])])]

testEq :: forall a.Reader (ValSet Double) (Equation (Position -> [a])) -> Equation a
testEq eq = getDiscEqInstance ( runReader eq $! initialGrid) testPosition


--writeTermsOrig:: (Num a, Show a, Fractional a)=> [Term a] -> String
writeTermsOrig terms vs=
    let writeTerm prev t= prev ++ case t of
            Unknown u -> show u ++ "X + "
            Constant c -> show c ++ " + "
            SubExpression s -> writeTerms (getTerms s) vs ++ " + "
            Derivative_adj d _ side _->
                "approxDeriv ( " ++ writeTerms [approximateDerivative_adj t testPosition vs] vs
                    ++ show d ++ " " ++ show  side ++" ) + "
    in foldl' writeTerm " " (reverse terms )

writeTerms :: [Term Double] -> AdjGraph -> String
writeTerms terms vs=
    let (_:_:xs) =  reverse $! writeTermsOrig terms vs
    in reverse xs

testPosition =   Position [1, 1, 0] 3 0
 -}

makeRows :: [[Double]] -> [Double]-> [Double] -> Int -> Int-> [[Double]]
makeRows whole curr [] _ _ = whole ++ [curr]
makeRows whole curr items 0 width = makeRows (whole ++ [curr] ) [] items width width
makeRows whole curr (x:xs) r width= makeRows whole (curr++[x]) xs (r-1) width

adjGraphToGrid vs prp rowLength zLevel=
    makeRows [] []
        (map (\(c,_) -> prop_adj Nondirectional prp c Prev vs )
            $ filter (\(_,b)-> getPositionComponent (origPos b) Z == zLevel )
                $ map (\x-> (x,getNode vs x) ) $ allPositionsCurrentTime vs
        )
        rowLength
        rowLength

--stringDomain:: (Num a, Fractional a, Show a ) => Property ->[Position]->Int-> ValSet a -> String
stringDomain prp  rowLength set zLevel =
    let rows = adjGraphToGrid set prp rowLength zLevel
        strRows = map ( foldl' (\prev next-> prev ++ " " ++ show next) "" ) rows
    in foldl' (\prev next -> prev ++ "\n" ++ next ) "" strRows

main:: IO AdjGraph
main =
    putStrLn "starting ..... "
    {-
    >>= (\_-> print ( solveUnknown testEquation ( Position [0, 0, 0] 3 0) initialGrid ))
    >>= (\_ -> putStrLn $ writeTerms (distributeMultiply testTerms 2) initialGrid)
    >>= (\_ -> print $ prop Nondirectional U testPosition Center  initialGrid)
    >>= (\_ -> putStrLn " continuity ------------ ")
    >>= (\_ -> putStrLn $ writeTerms ( rhs $ testEq continuity ) initialGrid )
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms ( lhs $ testEq continuity) initialGrid )
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> print $ solveUnknown (testEq continuity) testPosition initialGrid )
    >>= (\_ -> putStrLn " U Momentum------------ ")
    >>= (\_ -> putStrLn $ writeTerms ( rhs $ testEq uMomentum ) initialGrid )
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms ( lhs $ testEq uMomentum ) initialGrid )
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> print $ solveUnknown (testEq uMomentum) testPosition initialGrid )
    >>= (\_ -> putStrLn " V Momentum------------ ")
    >>= (\_ -> putStrLn $ writeTerms ( rhs $ testEq vMomentum ) initialGrid )
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms ( lhs $ testEq vMomentum ) initialGrid )
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> print $ solveUnknown (testEq vMomentum) testPosition initialGrid )
    >>= (\_ -> putStrLn " W Momentum------------ ")
    >>= (\_ -> putStrLn $ writeTerms ( rhs $ testEq wMomentum ) initialGrid )
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms ( lhs $ testEq wMomentum ) initialGrid )
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> print $ solveUnknown (testEq wMomentum) testPosition initialGrid )
    >>= (\_ -> putStrLn " ENERGY ------------ ")
    >>= (\_ -> putStrLn $ writeTerms ( rhs $ testEq energy ) initialGrid )
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms ( lhs $ testEq energy ) initialGrid )
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> print $ solveUnknown (testEq energy) testPosition initialGrid )
    -- >>= (\_ -> putStrLn $! stringDomain Pressure ( timePos $ offsetPosition (head $ calculatedPositions runTimeSteps) Prev) (1+maxPos X) runTimeSteps  )
    -- >>= (\_ -> putStrLn $! stringDomain Pressure (timePos $ head $ calculatedPositions runTimeSteps) (1+maxPos X) runTimeSteps  )
    -- >>= (\_ -> putStrLn $! stringDomain U ( timePos $ offsetPosition (head $ calculatedPositions runTimeSteps) Prev) (1+maxPos X) runTimeSteps  )
    -- >>= (\_ -> putStrLn $! stringDomain U (timePos $ head $ calculatedPositions runTimeSteps) (1+maxPos X) runTimeSteps  )
    -}
    -- >>= (\_ -> putStrLn $ show $  applyDiffEq_adj calcSteps_adj True initialGrid_adj)
    -- >>= (\_ -> print (totalDensity $ runSingleStep $ runSingleStep initialGrid_adj))
     >>= (\_ -> runTimeSteps_Print )

