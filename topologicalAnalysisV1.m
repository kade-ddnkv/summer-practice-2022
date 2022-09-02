% Предположим, мне нужно сравнить Pre_Relaxed и Post_Relaxed

% Загружаю уже предобработанные в EEGLAB данные.
% Эти данные уже переведены в структуру данных для fieldTrip.
savedir = ["C:\MATLAB\EEGData\MySaveFieldTrip"];
cd(savedir);
load('ft_EEG_RS_EO.mat')

EEG = ft_EEG.ASMRR.Run1;
clear ft_EEG;

% Каждый канал для ТАД я рассматриваю как точку в пространстве R^N, где N - это количество измерений на этом канале.
% То есть ломаную ЭЭГ на этом канале я рассматриваю как одну точку в R^N.
% Расстояние между каналами я определяю как стандартизированное (нормализованное) евклидово расстояние.

% Вычисляю стандартные отклонения перед основным кодом для ускорения.

pptnames = fieldnames(EEG)

for i_ppt = 1:length(pptnames)
	pptname = pptnames{i_ppt}
	EEG.(pptname).std = cell(1, height(EEG.(pptname).trialinfo));
	for i_trial = 1:height(EEG.(pptname).trialinfo)
		EEG.(pptname).std{i_trial} = std(EEG.(pptname).trial{1,i_trial});
	end
end

statedescript = {'Pre_Baseline', 'Post_BaseLine', 'Pre_Relaxed', 'Post_Relaxed'};

% Выбирается состояние.
state = 'Pre_Relaxed'

% Осторожно! Вычислене расстояний может работать долго.
% Для Pre_Relaxed код работал около 6 часов.
% Для Post_Relaxed код работал около 12 часов.
% Как идея для оптимизации: распараллелить вычисления (главное тут - учесть все проблемы многопоточного кода).

for i_ppt = 1:length(pptnames)
	pptname = pptnames{i_ppt}
	norm_euc_dists.(pptname).(state) = zeros(64,64);
	for band_from = 1:64
		for band_to = 1:(band_from - 1)
			temp1 = zeros(1,height(EEG.(pptname).trialinfo));
			for i_trial = 1:height(EEG.(pptname).trialinfo)
				triallabel = EEG.(pptname).trialinfo{i_trial,4}{1};
				if strcmp(triallabel, state) == 1
					temp2 = zeros(1, width(EEG.(pptname).trial{1,i_trial}));
					for k = 1:width(EEG.(pptname).trial{1,i_trial})
						temp2(1,k) = ((EEG.(pptname).trial{1,i_trial}(band_from,k) - EEG.(pptname).trial{1,i_trial}(band_to,k)) / EEG.(pptname).std{i_trial}(k));
					end
					temp1(1,i_trial) = sum(temp2.^2);
				end
			end
			norm_euc_dists.(pptname).(state)(band_from, band_to) = sqrt(sum(temp1));
		end
	end
	% Отзеркаливание относительно диагонали
	for y = 1:64
		for x = (y + 1):64
			norm_euc_dists.(pptname).(state)(y,x) = norm_euc_dists.(pptname).(state)(x,y);
		end
	end
end

% Можно строить картинки расстояний между каналами.
%cd img_norm_euc_dists
%pptname = 'P081'
%state = 'Post_Relaxed'
%imagesc(norm_euc_dists.(pptname).(state))
%set(gcf,'position',[10,10,600,510])
%title([pptname '-' strrep(state, '_', '-')])
%xlabel('Channel ID')
%ylabel('Channel ID')
%colorbar
%saveas(gcf,['norm_euc_dists_' pptname '_' state '.png'])
%cd ..

% Также выполняю предыдущий кусок кода и для Post_Relaxed.

save('norm_euc_dists', 'norm_euc_dists')
clearvars -except EEG norm_euc_dists pptnames

cd ('C:\Users\cyril\Desktop\Домашние работы\Практика 2 курс\Материалы для практики\My analysis code\Topological analysis\matlab_examples')
load_javaplex
cd ('..')


% Создание топологии.


max_dimension = 3;
% max_filtration_value = 2000; - это примерный предел для моего компьютера. В среднем вычисление одного человека занимает 3 минуты. Значит для Pre_Relaxed + Post_Relaxed будет около 2 часов.
max_filtration_value = 2000;
num_divisions = 1000;


pptnames = fieldnames(EEG);

states = {'Pre_Relaxed', 'Post_Relaxed'}

% Фильтрация Виеториса-Рипса.

for i_state = 1:length(states)
	state = states{i_state}
	for i_ppt = 1:length(pptnames)
		pptname = pptnames{i_ppt}
		
		point_cloud = norm_euc_dists.(pptname).(state);

		% create a Vietoris-Rips stream 
		stream = api.Plex4.createVietorisRipsStream(point_cloud, max_dimension, max_filtration_value, num_divisions);
		num_simplices = stream.getSize()

		% get persistence algorithm over Z/2Z
		persistence = api.Plex4.getModularSimplicialAlgorithm(max_dimension, 2);

		% compute the intervals
		intervals.(state).(pptname) = persistence.computeAnnotatedIntervals(stream);
	end
end
% intervals не могут быть сериализованы, поэтому не сохраняются.


% Мы получили штрихкоды персистентности
% Из этих штрихкодов нужно достать те каналы (вершины), которые образуют дырки в топологии.
% Документация по интервалам в javaPlex:
% https://appliedtopology.github.io/javaplex/doc/allclasses-noframe.html
% https://github.com/appliedtopology/javaplex/issues/11
% https://github.com/appliedtopology/primitive-lib/blob/master/src/edu/stanford/math/primitivelib/autogen/formal_sum/IntSparseFormalSum.java

% Здесь я из штрихкодов достаю:
% 1) Длины интервалов
% 2) Частоту того, сколько раз каждый канал участвовал в образовании дыр в топологии.

for i_state = 1:length(states)
	state = states{i_state}
	interval_lens.(state) = [];
	vertices_freq.(state).full = java.util.HashMap;
	for i_ppt = 1:length(pptnames)
		pptname = pptnames{i_ppt}
		vertices_freq.(state).(pptname) = java.util.HashMap;
		
		% Нахожу точки, которые образуют в топологии пустоты размерности dim.
		% Размерность 0 не интересует.
		for dim = 1:(max_dimension-1)
			gens = intervals.(state).(pptname).getIntervalGeneratorPairsAtDimension(dim);
			for i_gen = 0:(gens.size-1)
				interval = gens.get(i_gen).getFirst();
				if (interval.isRightInfinite())
					interval_len = max_filtration_value - cell2mat(cell(interval.getStart()));
				else
					interval_len = cell2mat(cell(interval.getEnd())) - cell2mat(cell(interval.getStart()));
				end
				interval_lens.(state) = [interval_lens.(state) interval_len];
				
				it = gens.get(i_gen).getSecond().iterator();
				gen_v_set = java.util.HashMap;
				while (it.hasNext())
					it.advance();
					v = it.key().getVertices();
					for i_v = 1:length(v)
						gen_v_set.put(v(i_v), 1);
					end
				end
				gen_v_set = gen_v_set.keySet().toArray();
				for i_v = 1:length(gen_v_set)
					v = gen_v_set(i_v);
					if (vertices_freq.(state).(pptname).containsKey(v))
						vertices_freq.(state).(pptname).put(v, vertices_freq.(state).(pptname).get(v) + interval_len);
					else
						vertices_freq.(state).(pptname).put(v, interval_len);
					end
					if (vertices_freq.(state).full.containsKey(v))
						vertices_freq.(state).full.put(v, vertices_freq.(state).full.get(v) + interval_len);
					else
						vertices_freq.(state).full.put(v, interval_len);
					end
				end
			end
		end
	end
end

save('vertices_freq', 'vertices_freq')
clearvars -except EEG norm_euc_dists pptnames states intervals interval_lens vertices_freq max_dimension max_filtration_value num_divisions


% ПЛОХОЙ ПРОШЛЫЙ КОД

% Строю гистограммы (Для каждого канала указывается частота, с которой он участвует в создании пустот размерности >=1 в топологии).

%states = fieldnames(vertices_freq);

%for i_state = 1:length(states)
%	state = states{i_state}
%	channels = cell2mat(cell(vertices_freq.(state).keySet.toArray));
%	[channels_sorted, channels_order] = sort(channels);
%	freqs = cell2mat(cell(vertices_freq.(state).values.toArray));
%	freqs = freqs(channels_order);
%	bar(freqs)
%	title(strrep(state, '_', '-'))
%	xlabel('Channel ID')
%	ylabel('Frequency')
%	saveas(gcf,['chan_freq_' state '.png'])
%end


% Картинки Pre_Relaxed и Post_Relaxed одинаковые.

% Где я мог ошибиться:
% 1) Я считал, что выборки Pre_Relaxed и Post_Relaxed являются независимыми. Однако это не так. Каждый человек (почти каждый) имеет триалы и с Pre_Relaxed, и с Post_Relaxed, то есть используется в двух выборках. Не уверен, могло ли это повлиять.
% 2) Я брал любые дырки в топологии, и существующие долго, и почти сразу же умирающие. Это, наверное, неправильно. То есть нужно понять, как отделить значимые интервалы от незначимых.


% Я понял, что я делал неправильно. Я учитывал то, что интервал (жизни дырки) образуют вершины, но не учитывал то, сколько этот интервал длится. А это важно, ведь чем больше время жизни дырки, тем важнее она.
% Другими словами, раньше я для каждого канала вычислял частоту в кол-ве дырок, а нужно вычислять в кол-ве тактов. Такт - это одна итерация в фильтрации Виеториса-Рипса. В этом случае я должен вес каждой вершины просто умножать на длину интервала дырки, в создании которой она участвовала.


% Для начала я сделаю гистограмму длин интервалов, на всякий случай, посмотреть.

cd img_interval_lengths
for i_state = 1:length(states)
	state = states{i_state}
	h = histogram(interval_lens.(state), 70)
	title(strrep(state, '_', '-'))
	xlabel('Interval Length')
	ylabel('Frequency')
	saveas(gcf,['histogram_interval_lengths_' state '.png'])
end
cd ..

% Можно раскомментировать, чтобы посмотреть на баркоды и числа Бетти.

%state = 'Post_Relaxed';
%pptname = 'P081';
%cd img_barcodes
%% create the barcode plots
%options.filename = [pptname '-' strrep(state, '_', '-')];
%options.max_filtration_value = max_filtration_value;
%options.max_dimension = max_dimension - 1;
%options.side_by_side = true;
%plot_barcodes(intervals.(state).(pptname), options);
%cd ..

%% get the infinite barcodes
%infinite_barcodes = intervals.getInfiniteIntervals();

%% print out betti numbers array
%betti_numbers_array = infinite_barcodes.getBettiSequence()


% Что я пытаюсь сделать:
% Определить, какие каналы как участвуют в создании дырок.
% Например, может оказаться, что в экспериментальных условиях (Post_Relaxed) некоторые каналы гораздо чаще создают длительные дырки в топологии, чем в базовых условиях (Pre_Relaxed).

cd img_chan_freq
states = fieldnames(vertices_freq);
for i_state = 1:length(states)
	state = states{i_state}
	channels = cell2mat(cell(vertices_freq.(state).keySet.toArray));
	[channels_sorted, channels_order] = sort(channels);
	freqs = cell2mat(cell(vertices_freq.(state).values.toArray));
	freqs = freqs(channels_order);
	bar(freqs)
	title(strrep(state, '_', '-'))
	xlabel('Channel ID')
	ylabel('Frequency (in Rips ticks)')
	saveas(gcf,['chan_freq_' state '.png'])
end
cd ..


% Итого:
% Пики на гистограммах в одинаковых местах по горизонтали.
% Но гистограмма Post_Relaxed как будто выросла раза в 1.5-2 по сравнению с гистограммой Pre_Relaxed.
% Это не может быть обусловлено выбросами (одним-двумя людьми, у которых на канале такое случилось): можно посмотреть на гистограммы длин интервалов - там видно, что максимальная длина интервала у Post_Relaxed - это 600, и то это было у одного человека в один момент.

% Возникает такое чувство, что в Post_Relaxed все дырки в топологии держатся в 2 раза дольше.




% Как понять, эти различия существенны или нет?
% Я недавно узнал про пермутационный тест и считаю, что он тут хорошо подойдет.

% Можно вычислить такие гистограммы для каждого участника и сравнивать выборки по каждому каналу. Типо, в Pre_Relaxed для первого канала у нас есть примерно 16 частот (всего 25 участников), они образуют выборку, и также есть выборка для Post_Relaxed. Проводим пермутационный тест, узнаем, значимо ли отличие в этом канале или нет.

% Так и сделаю.


significant_chans = [];
for i_chan = 1:64
	for i_state = 1:length(states)
		state = states{i_state};
		% Выборка из 25 элементов.
		selection.(state) = [];
		for i_ppt = 1:length(pptnames)
			pptname = pptnames{i_ppt};
			selection.(state) = [selection.(state) vertices_freq.(state).(pptname).get(i_chan - 1)];
		end
	end
	% Можно добавить 'plotresult', 1 в конец функции, чтобы выводило гистограмму различий средних.
	[p observeddifference] = permutationTest(selection.Pre_Relaxed, selection.Post_Relaxed, 10000);
	chan.index = i_chan;
	chan.p = p;
	chan.observeddifference = observeddifference;
	if (p <= 0.05)
		significant_chans = [significant_chans chan];
	end
end

length(significant_chans) % 52
% Для 52 каналов различия оказались значимыми.
% То есть почти для всех.

% Гистограмма отличий средних:
histogram(arrayfun(@(x) x.observeddifference, significant_chans), 20)
saveas(gcf,['observeddifference_significant_chans.png'])
% В среднем каналы на 140 тактов (инкрементов фильтрации Виеториса-Рипса) дольше держат дырки в топологии.

% Вывод.
% На диаграмме частот участия каналов в образовании дыр видно, что пики находятся в одних и тех же каналах.
% Но в условиях Post_Relaxed у людей дырки в топологии держатся дольше.
% Помню, что дистанция между каналами рассчитывается как нормированное евклидово расстояние.
